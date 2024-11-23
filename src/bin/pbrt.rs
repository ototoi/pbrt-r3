use clap::*;

use pbrt_r3::core::api::*;
use pbrt_r3::core::parser::*;
use pbrt_r3::core::pbrt::*;
use pbrt_r3::core::stats;
use pbrt_r3::displays::TevDisplay;
use std::cell::RefCell;
use std::env;
use std::path::Path;
use std::path::PathBuf;
use std::process;
use std::sync::Arc;
use std::sync::RwLock;
use std::thread::available_parallelism;

use log::*;
/*
  --cropwindow <x0,x1,y0,y1> Specify an image crop window.
  --help               Print this help text.
  --nthreads <num>     Use specified number of threads for rendering.
  --outfile <filename> Write the final image to the given filename.
  --quick              Automatically reduce a number of quality settings to
                       render more quickly.
  --quiet              Suppress all text output other than error messages.

Logging options:
  --logdir <dir>       Specify directory that log files should be written to.
                       Default: system temp directory (e.g. $TMPDIR or /tmp).
  --logtostderr        Print all logging messages to stderr.
  --minloglevel <num>  Log messages at or above this level (0 -> INFO,
                       1 -> WARNING, 2 -> ERROR, 3-> FATAL). Default: 0.
  --v <verbosity>      Set VLOG verbosity.

Reformatting options:
  --cat                Print a reformatted version of the input file(s) to
                       standard output. Does not render an image.
  --toply              Print a reformatted version of the input file(s) to
                       standard output and convert all triangle meshes to
                       PLY files. Does not render an image.
*/

#[derive(Debug, Parser)]
#[clap(author, about, version, disable_help_flag = true)]
struct CommandOptions {
    /// Input .pbrt file.
    #[arg(short, long, value_name = "filename")]
    pub infile: Option<PathBuf>,

    /// Write the final image to the given filename.
    #[arg(short, long, value_name = "filename")]
    pub outfile: Option<PathBuf>,

    /// Specify an image crop window.
    #[arg(long, value_delimiter = ',', value_name = "x0,x1,y0,y1")]
    pub cropwindow: Option<Vec<f32>>,

    /// Print this help text.
    #[arg(short, long, action = clap::ArgAction::HelpLong)]
    pub help: Option<bool>,

    /// Use specified number of threads for rendering.
    #[arg(short = 'j', long = "nthreads", value_name = "num")]
    pub nthreads: Option<i32>,

    /// Automatically reduce a number of quality settings to render more quickly.
    #[arg(long, default_value = "false")]
    pub quick: bool,

    /// Suppress all text output other than error messages.
    #[clap(long, default_value = "false")]
    pub quiet: bool,

    // Logging options
    /// Specify directory that log files should be written to.
    /// Default: system temp directory (e.g. $TMPDIR or /tmp).
    #[arg(long, value_name = "dir")]
    pub logdir: Option<PathBuf>,

    /// Print all logging messages to stderr.
    #[arg(long, default_value = "false")]
    pub logtostderr: bool,

    /// Log messages at or above this level (0 -> INFO,
    /// 1 -> WARNING, 2 -> ERROR, 3-> FATAL).
    #[arg(long, value_name = "num")] //value_enum
    pub minloglevel: Option<i32>,

    /// Set VLOG verbosity.
    //#[arg(long = "v", default_value = "0", value_name = "verbosity")]
    //pub verbosity: i32,

    // Reformatting options
    /// Print a reformatted version of the input file(s) to standard output.
    /// Does not render an image.
    #[arg(short, long, default_value = "false")]
    pub cat: bool,

    /// Print a reformatted version of the input file(s) to standard output and convert all triangle meshes to PLY files.
    /// Does not render an image.
    #[arg(short, long, default_value = "false")]
    pub toply: bool,

    /// Display-server ex. localhost:14158
    #[arg(long = "display-server", value_name = "url")]
    pub display_server: Option<String>,

    /// Set Pixelsamples.
    #[arg(short = 's', long = "pixelsamples", value_name = "num")]
    pub pixelsamples: Option<i32>,

    #[arg(value_name = "filename.pbrt")]
    pub pbrtfile: Option<Vec<PathBuf>>,
}

fn init_logger(opts: &CommandOptions) {
    if let Some(minloglevel) = opts.minloglevel {
        const LOG_LEVELS: &[&str] = &["trace", "debug", "info", "warn", "error"];
        let log_level = LOG_LEVELS[(minloglevel + 2).clamp(0, 4) as usize];
        env::set_var("RUST_LOG", log_level);
    } else {
        //default log level : warn
        let log_level = env::var("RUST_LOG").unwrap_or_else(|_| "warn".to_owned());
        env::set_var("RUST_LOG", log_level);
    }

    env_logger::Builder::from_default_env()
        //.format_timestamp(None)
        .format_target(false)
        .format_module_path(false)
        .init();
}

fn print_scene(input_path: &Path, _opts: &CommandOptions) -> i32 {
    let mut context = MutipleContext::new();
    context.add(Arc::new(RefCell::new(PrintContext::new_stdout(false))));

    let input_path = String::from(input_path.to_str().unwrap());
    //should
    if pbrt_parse_file(&input_path, &mut context).is_ok() {
        0
    } else {
        -1
    }
}

fn create_scene_and_integrator(
    input_path: &Path,
    opts: &CommandOptions,
) -> Result<(Arc<Scene>, Arc<RwLock<dyn Integrator>>), PbrtError> {
    let mut context = MutipleContext::new();

    let scene_context = Arc::new(RefCell::new(SceneContext::new()));
    context.add(scene_context.clone());

    let input_path = String::from(input_path.to_str().unwrap());

    pbrt_parse_file(&input_path, &mut context)?;

    if let Some(pixelsamples) = opts.pixelsamples {
        let pixelsamples = i32::max(1, pixelsamples);
        let mut sc = scene_context.borrow_mut();
        sc.replace_params_int("Sampler", "pixelsamples", pixelsamples)?;
    }

    if let Some(outfile) = opts.outfile.as_ref() {
        let outfile = String::from(outfile.to_str().unwrap());
        let mut sc = scene_context.borrow_mut();
        sc.replace_params_string("Film", "filename", &outfile)?;
    }

    let sc = scene_context.borrow();
    let scene = sc.make_scene();
    let integrator = sc.make_integrator();
    if let Some(s) = scene {
        if let Some(i) = integrator {
            return Ok((s, i));
        }
    }
    Err(PbrtError::error("Coudn't read file."))
}

fn create_display(hostname: &str) -> Result<Arc<RwLock<dyn Display>>, PbrtError> {
    let mut tev = TevDisplay::new();
    tev.connect(hostname)?;
    Ok(Arc::new(RwLock::new(tev)))
}

fn render_scene(input_path: &Path, opts: &CommandOptions) -> i32 {
    stats::clear_stats();
    if !opts.quiet {
        let nthreads = available_parallelism().unwrap().get();
        let version = env!("CARGO_PKG_VERSION");
        println!("pbrt-r3 version {} [Detected {} cores]", version, nthreads);
        println!();
        println!("This is an unofficial Rust port of pbrt-v3.");
        println!("The source code of this port is covered by the BSD License.");
        println!("See the file LICENSE.txt for the conditions of the license.");
        println!();
        println!("The license for the original implementation pbrt-v3 is as follows:");
        println!("--------------------------------------------------------------------------------------");
        println!("Copyright (c)1998-2018 Matt Pharr, Greg Humphreys, and Wenzel Jakob.");
        println!(
            "The source code to pbrt (but *not* the book contents) is covered by the BSD License."
        );
        println!("See the file LICENSE.txt for the conditions of the license.");
        println!("--------------------------------------------------------------------------------------");
        println!();
    }
    let r = create_scene_and_integrator(input_path, opts);
    match r {
        Ok((scene, integrator)) => {
            if let Some(hostname) = opts.display_server.as_ref() {
                match create_display(hostname) {
                    Ok(display) => {
                        let integrator = integrator.as_ref().read().unwrap();
                        let camera = integrator.get_camera();
                        let film = camera.as_ref().get_film();
                        let mut f = film.as_ref().write().unwrap();
                        f.add_display(&display);
                    }
                    Err(e) => {
                        warn!("{}", e);
                    }
                }
            }

            if !opts.quiet {
                stats::print_stats();
            }

            {
                let mut integrator = integrator.as_ref().write().unwrap();
                let scene = scene.as_ref();
                integrator.render(scene);
            }
        }
        Err(e) => {
            let msg = format!("{}", e);
            error!("{}", msg);
            return -1;
        }
    }

    println!("\n");
    if !opts.quiet {
        stats::print_stats();
        stats::clear_stats();
    }

    return 0;
}

pub fn main() {
    let opts = CommandOptions::parse();
    init_logger(&opts);
    let input = if let Some(infiles) = opts.pbrtfile.as_ref() {
        Some(infiles[0].clone())
    } else {
        opts.infile.as_ref().cloned()
    };

    if let Some(ipath) = input.as_ref() {
        if !ipath.exists() {
            println!("{}", CommandOptions::command().render_usage());
            process::exit(-1);
        }
    } else {
        println!("{}", CommandOptions::command().render_usage());
        process::exit(-1);
    }

    let input_path = input.unwrap();
    let input_path = input_path.canonicalize().unwrap();

    if opts.cat {
        let ret = print_scene(&input_path, &opts);
        process::exit(ret);
    } else if opts.toply {
        //todo: implement toply
    } else {
        let ret = render_scene(&input_path, &opts);
        process::exit(ret);
    }
}

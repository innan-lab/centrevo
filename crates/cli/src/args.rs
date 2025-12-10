use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct InitArgs {
    /// Simulation name
    #[arg(short = 'N', long, default_value = "simulation")]
    pub name: String,

    /// Output database path
    #[arg(short, long, default_value = "simulation.db")]
    pub output: PathBuf,

    /// Population size
    #[arg(short = 'n', long, default_value = "100")]
    pub population_size: usize,

    /// Number of generations
    #[arg(short = 'g', long, default_value = "1000")]
    pub generations: usize,

    /// Repeat unit length
    #[arg(long, default_value = "171")]
    pub ru_length: usize,

    /// Repeat units per HOR
    #[arg(long, default_value = "12")]
    pub rus_per_hor: usize,

    /// HORs per chromosome
    #[arg(long, default_value = "100")]
    pub hors_per_chr: usize,

    /// Chromosomes per haplotype
    #[arg(long, default_value = "1")]
    pub chrs_per_hap: usize,

    /// Mutation rate (Uniform/Jukes-Cantor)
    ///
    /// Defaults to 1e-5. Mutually exclusive with specific rates.
    #[arg(
        long,
        conflicts_with_all = ["rate_ac", "rate_ag", "rate_at", "rate_cg", "rate_ct", "rate_gt"]
    )]
    pub mutation_rate: Option<f64>,

    /// Mutation rate A<->C
    #[arg(long)]
    pub rate_ac: Option<f64>,
    /// Mutation rate A<->G
    #[arg(long)]
    pub rate_ag: Option<f64>,
    /// Mutation rate A<->T
    #[arg(long)]
    pub rate_at: Option<f64>,
    /// Mutation rate C<->G
    #[arg(long)]
    pub rate_cg: Option<f64>,
    /// Mutation rate C<->T
    #[arg(long)]
    pub rate_ct: Option<f64>,
    /// Mutation rate G<->T
    #[arg(long)]
    pub rate_gt: Option<f64>,

    /// Insertion rate (per base per generation)
    #[arg(long, default_value = "0.0")]
    pub indel_ins_rate: f64,

    /// Deletion rate (per base per generation)
    #[arg(long, default_value = "0.0")]
    pub indel_del_rate: f64,

    /// Geometric parameter p for indel length (higher p = shorter indels)
    #[arg(long, default_value = "0.5")]
    pub indel_length_p: f64,

    /// Recombination break probability
    ///
    /// Defaults to 1e-6 or approx. 1 DSB per Mb per generation.
    #[arg(long, default_value = "1e-6")]
    pub recomb_rate: f64,

    /// Crossover probability (given break)
    ///
    /// Defaults to 0.01. Crossovers are heavily suppressed in centromeres; most events are gene conversions.
    #[arg(long, default_value = "0.01")]
    pub crossover_prob: f64,

    /// Gene conversion extension probability
    ///
    /// Defaults to 0.95. Results in average tract length ~20bp (1/(1-0.95)).
    #[arg(long, default_value = "0.95")]
    pub gc_extension_prob: f64,

    /// Homology strength (0.0 = random, >0.0 = preference for similarity)
    ///
    /// Defaults to 5.0. (Strong preference for homologous sequences to drive homogenization)
    #[arg(long, default_value = "5.0")]
    pub homology_strength: f64,

    /// Search window for homology (in RUs)
    ///
    /// Defaults to 100. (Allows interaction with neighboring ~10-20kb of sequence)
    #[arg(long, default_value = "100")]
    pub search_window: usize,

    /// Optimal GC content for fitness (0.0-1.0)
    #[arg(long, requires = "fit_gc_conc")]
    pub fit_gc_opt: Option<f64>,

    /// Concentration parameter for GC fitness
    #[arg(long, requires = "fit_gc_opt")]
    pub fit_gc_conc: Option<f64>,

    /// Optimal chromosome length for fitness (bases)
    #[arg(long, requires = "fit_len_std")]
    pub fit_len_opt: Option<usize>,

    /// Standard deviation for length fitness (bases)
    #[arg(long, requires = "fit_len_opt")]
    pub fit_len_std: Option<f64>,

    /// Shape parameter for sequence similarity fitness
    #[arg(long)]
    pub fit_seq_sim: Option<f64>,

    /// Shape parameter for length similarity fitness
    #[arg(long)]
    pub fit_len_sim: Option<f64>,

    /// Default recording interval (record every N generations)
    #[arg(long, default_value = "100")]
    pub record_every: usize,

    /// Random seed
    #[arg(long)]
    pub seed: Option<u64>,
}

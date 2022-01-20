{ bionix, ref, binWidth ? 1, blacklist ? null, targets }:

with bionix;
with lib;
with pkgs;

let
  r = rWrapper.override {
    packages = with rPackages;[ QDNAseq Biobase (callPackage ./bsgenome.nix { inherit ref; }) ];
  };

  sequences = map (x: head (splitString ":" x)) targets;

  runRScript = name: script: stage {
    inherit name;
    buildInputs = [ r ];
    buildCommand = ''
      Rscript ${writeText "${name}.R" script}
    '';
  };

in
makeExtensible (self: with self;
{
  toBW = exec
    (_: wig: stage {
      name = "wig2bigwig.bw";
      buildInputs = [ pkgs.kent ];
      buildCommand = ''
        wigToBigWig ${wig}/*.{wig,chrom.sizes} $out
      '';
    });

  bins = exec'' (runRScript "bin" ''
    library(BSgenome.vivax)
    library(QDNAseq)
    library(Biobase)
    pfBins <- createBins(BSgenome.vivax, ${toString binWidth})
    pfBins$mappability <- calculateMappability(pfBins
      , bigWigAverageOverBed="${kent}/bin/bigWigAverageOverBed"
      , bigWigFile="${toBW {} (bionix.genmap.calcmap {} ref)}"
      , chrPrefix="")
    ${optionalString (blacklist != null) ''
      pfBins$blacklist <- calculateBlacklist(pfBins, bedFiles="${blacklist}")
    ''}
    saveRDS(pfBins, Sys.getenv("out"))
  '');

  count = flip def { mem = 10; walltime = "12:00:00"; } (exec ({ name ? null, ... }: input: runRScript "count" ''
    library(QDNAseq)
    pfBins <- readRDS("${bins}")
    seqnames <- c(${concatMapStringsSep "," (x: "'${x}'") sequences})
    pfBins2chrom <- pfBins[pfBins$chromosome %in% seqnames,]
    countsbin <- binReadCounts(pfBins2chrom, bamfiles="${input}" ${optionalString (name != null) ", bamnames='${name}'"})
    saveRDS(countsbin, Sys.getenv("out"))
  ''));

  call = flip def { mem = 10; walltime = "12:00:00"; } (exec (attrs: input: runRScript "call" ''
    library(QDNAseq)
    countsbin <- readRDS("${count attrs input}")
    countsFiltered <- applyFilters(
        countsbin
      , mappability=50
      , blacklist=${if blacklist == null then "FALSE" else "TRUE"}
      , residual=F)
    countsCorrected <- estimateCorrection(countsFiltered)
    copyNums <- correctBins(countsCorrected)
    copyNumsNormed <- normalizeBins(copyNums)
    saveRDS(copyNumsNormed, Sys.getenv("out"))
  ''));
})

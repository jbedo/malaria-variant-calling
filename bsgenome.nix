{ rWrapper, rPackages, ref, writeText, runCommand }:

let

  r = rWrapper.override {
    packages = with rPackages; [ BSgenome ];
  };

  forgeBSgenome = ref:
    let
      script = writeText "forgeBSgenome.r" ''
        library("BSgenome")
        forgeBSgenomeDataPkg("${seed}", "${splitRef}", getwd())
      '';
    in
    runCommand "forgeBSgenome" { nativeBuildInputs = [ r ]; } ''
      Rscript ${script}
      cp -r BSgenome.vivax $out
    '';

  seed = writeText "seed.txt" ''
    Package: BSgenome.vivax
    Title: Plasmodium vivax
    Description: Plasmodium vivax
    Version: 0.1
    Author: John Smith
    Maintainer: John Smith <john@smith.com>
    License: public domain
    organism: Plasmodium vivax P01
    common_name: P. vivax
    provider: PlasmoDB
    genome: PvivaxP01
    release_date: 1970-01-01
    organism_biocview: Plasmodium_vivax_P01
    BSgenomeObjname: PvivaxP01
    source_url: ${builtins.head ref.urls}
    circ_seqs: character(0)
    seqs_srcdir: ${splitRef}/
    seqfiles_suffix:
    seqnames: dir("${splitRef}")
  '';

  splitRef = runCommand "splitFA" { } ''
    mkdir $out
    cd $out
    awk '/^>/{output=substr($1,2,length($1)-1)}{print > output}' ${ref}
  '';

in

rPackages.buildRPackage {
  name = "BSGenome.Vivax";
  src = forgeBSgenome ref;
  propagatedBuildInputs = with rPackages; [ BSgenome ];
}

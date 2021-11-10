{ bionix ? import <bionix> { } }:

with bionix;
with lib;
with types;

let
  ref = fetchFastA {
    url = "https://plasmodb.org/common/downloads/Current_Release/PvivaxP01/fasta/data/PlasmoDB-54_PvivaxP01_Genome.fasta";
    sha256 = "sha256-1ZahI54i1iZ3kvJGzEz00Gv3l9BgzHL+1gJ1I+1wA3k=";
  };
  PMIX = "PvP01_13_v2:876959-881753";
  PMX = "PvP01_01_v2:555499-559411";

  ena = pkgs.callPackage ./ena.nix { inherit fetchFastQGZ; };

  samples = ena.importFromJSON ./PRJEB2140.json;

  pipeline = flip pipe [
    (mapAttrsToList (id: bwa.align { inherit ref; flags = "-R'@RG\\tID:${id}\\tSM:${id}'"; }))
    (map (samtools.sort { }))
    (octopus.call { targets = [ PMIX PMX ]; flags = "-P 1"; })
  ];

in
pipeline samples

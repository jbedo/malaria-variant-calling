self: super:

with super;
with lib;
with types;

let
  getRef = matchFiletype "gatk" { gz = getRef; bam = { ref, ... }: ref; vcf = { ref, ... }: ref; };
  lnRef = input: ''
    ln -s ${getRef input} ref.fa
    ln -s ${samtools.faidx {} (getRef input)} ref.fa.fai
    ln -s ${samtools.dict {} (getRef input)} ref.dict
  '';
  lnVcf = input: ''
    ln -s ${input} ./$(basename ${input})
    ln -s ${self.gatk.index {} input} ./$(basename ${input}).tbi
  '';
in
{
  gatk =
    {
      app = pkgs.callPackage ./gatk-app.nix { };

      callHaplotype = exec ({ targets, ... }: input: stage {
        name = "HaplotypeCaller.g.vcf.gz";
        buildInputs = [ self.gatk.app ];
        buildCommand = ''
          ${lnRef input}
          ln -s ${input} ./input.bam
          ln -s ${samtools.index {} input} ./input.bam.bai
          gatk HaplotypeCaller \
            -R ref.fa \
            -I input.bam \
            -L ${pkgs.writeText "targets.intervals" (concatStringsSep "\n" targets)} \
            -O $out \
            -ERC GVCF \
            -ploidy 1 \
            --native-pair-hmm-threads $NIX_BUILD_CORES
        '';
        passthru.filetype = filetype.vcf { ref = getRef input; };
      });

      merge = exec (_: xs: stage {
        name = "CombineGVCFs.g.vcf.gz";
        buildInputs = [ self.gatk.app ];
        buildCommand = ''
          ${lnRef (head xs)}
          ${concatMapStringsSep "\n" lnVcf xs}
          gatk CombineGVCFs \
            -R ref.fa \
            ${concatMapStringsSep "\n" (x: "--variant $(basename ${x}) \\") xs}
            -O $out
        '';
        passthru.filetype = filetype.vcf { ref = getRef (head xs); };
      });

      callGenotypes = exec (_: input: stage {
        name = "GenotypeGVCFs.vcf.gz";
        buildInputs = [ self.gatk.app ];
        buildCommand = ''
          ${lnRef input}
          ${lnVcf input}
          gatk GenotypeGVCFs \
            -R ref.fa \
            -V $(basename ${input})\
            -O $out
        '';
        passthru.filetype = filetype.vcf { ref = getRef input; };
      });

      index = exec
        (_: input: stage {
          name = "IndexFeatureFile.tbi";
          buildInputs = [ self.gatk.app ];
          buildCommand = ''
            gatk IndexFeatureFile \
             -I ${input} \
             -O $out
          '';
        });
    };
}

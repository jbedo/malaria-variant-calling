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
  lnVcf = name: input: ''
    ln -s ${input} ./${name}.vcf.gz
    ln -s ${self.gatk.index {} input} ./${name}.vcf.gz.tbi
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
            -ploidy 1
        '';
        passthru.filetype = filetype.vcf { ref = getRef input; };
        passthru.multicore = true;
      });

      merge = exec (_: { a, b }: assert getRef a == getRef b; stage {
        name = "CombineGVCFs.g.vcf.gz";
        buildInputs = [ self.gatk.app ];
        buildCommand = ''
          ${lnRef a}
          ${lnVcf "a" a}
          ${lnVcf "b" b}
          gatk CombineGVCFs \
            -R ref.fa \
            --variant a.vcf.gz \
            --variant b.vcf.gz \
            -O $out
        '';
        passthru.filetype = filetype.vcf { ref = getRef a; };
        passthru.multicore = true;
      });

      callGenotypes = exec (_: input: stage {
        name = "GenotypeGVCFs.vcf.gz";
        buildInputs = [ self.gatk.app ];
        buildCommand = ''
          ${lnRef input}
          ${lnVcf "input" input}
          gatk GenotypeGVCFs \
            -R ref.fa \
            -V input.vcf.gz \
            -O $out
        '';
        passthru.filetype = filetype.vcf { ref = getRef input; };
        passthru.multicore = true;
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

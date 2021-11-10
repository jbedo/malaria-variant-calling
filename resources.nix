self: super: with super;
{
  bwa.align = def bwa.align { mem = 20; ppn = 20; walltime = "1:00:00"; };
  samtools.sort = def samtools.sort { ppn = 10; mem = 15; flags = "-m 1G"; walltime = "1:00:00"; };
  octopus.call = def octopus.call { mem = 20; ppn = 24; walltime = "24:00:00"; };
}

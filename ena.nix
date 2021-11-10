{ lib, fetchFastQGZ }:

{
  importFromJSON = input: with lib; let
    meta = with builtins; pipe input [ readFile fromJSON ];
    fetch = { run_accession, fastq_md5, fastq_ftp, ... }:
      let
        fqs = splitString ";" fastq_ftp;
        md5s = splitString ";" fastq_md5;
        snd = x: head (tail x);
        fq1 = head fqs;
        fq2 = snd fqs;
        md51 = head md5s;
        md52 = snd md5s;
        go = url: hash: fetchFastQGZ {
          url = "ftp://" + url;
          outputHash = hash;
          outputHashAlgo = "md5";
        };
      in
      nameValuePair run_accession
        {
          input1 = go fq1 md51;
          input2 = go fq2 md52;
        };
  in
  listToAttrs (map fetch meta);
}

{ stdenvNoCC, fetchurl, jre, makeWrapper, unzip }:

stdenvNoCC.mkDerivation rec {
  pname = "gatk";
  version = "4.2.3.0";

  src = fetchurl {
    url = "https://github.com/broadinstitute/gatk/releases/download/4.2.3.0/gatk-4.2.3.0.zip";
    sha256 = "sha256-EvvQMUIxFBmgUNuoAJ1hXuj6zrMrRqrY9ESNr9YeEes=";
  };

  nativeBuildInputs = [ makeWrapper unzip ];

  buildPhase = false;

  installPhase = ''
    install -Dm 644 gatk-package-${version}-local.jar $out/libexec/gatk.jar
    makeWrapper ${jre}/bin/java $out/bin/gatk \
      --add-flags "-jar $out/libexec/gatk.jar"
  '';
}

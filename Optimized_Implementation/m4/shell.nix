{pkgs ? import <nixpkgs> {}, ...}: let
  envname = "pio-arduino-fhs";
  # as a function to make sure the same pkgs is used as in targetPkgs
  mypython = pks: pks.python3.withPackages (ps: with ps; [pks.platformio pylibftdi pyusb]);
  # "proxy" env, is this useful/necessary???
  # myEnv = pkgs.buildEnv {
  #   name = envname;
  #   paths = [ pkgs.zsh ];
  # };
in
  (pkgs.buildFHSUserEnv {
    name = envname;
    targetPkgs = pkgs: (with pkgs; [
      arduino-cli
      avrdude
      libftdi
      libftdi1
      libusb1
      platformio
      platformio-core
      (mypython pkgs)
      zsh
      vscode 
      clang 
      zstd
    ]);
  })
.env

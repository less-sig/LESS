with import <nixpkgs> {};
let 

  unstable_pkgs = import <nixos-unstable> {};
  bInputs = [ 
    # just for cargo 
    clang 
    llvm 
    lldb
    unstable_pkgs.cargo
    unstable_pkgs.rustc
  ] ++ (lib.optionals pkgs.stdenv.isLinux ([
  ]));
in
{ pkgs ? import <nixpkgs> {} }:

stdenv.mkDerivation {
  name = "";
  src = ./.;

  buildInputs = bInputs; #with pkgs; [ ];
  nativeBuildInputs = with pkgs; [];
  RUST_SRC_PATH = "${pkgs.rust.packages.stable.rustPlatform.rustLibSrc}";
}

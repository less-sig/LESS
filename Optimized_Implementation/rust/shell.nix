with import <nixpkgs> {};
let 
  bInputs = [ 
    # just for cargo 
    clang 
    llvm 
    lldb
    cargo
    rustc
  ] ++ (lib.optionals pkgs.stdenv.isLinux ([
  ]));
in
{ pkgs ? import <nixpkgs> {} }:

stdenv.mkDerivation {
  name = "";
  src = ./.;

  buildInputs = bInputs; #with pkgs; [ ];
  nativeBuildInputs = with pkgs; [addOpenGLRunpath];
  RUST_SRC_PATH = "${pkgs.rust.packages.stable.rustPlatform.rustLibSrc}";
}

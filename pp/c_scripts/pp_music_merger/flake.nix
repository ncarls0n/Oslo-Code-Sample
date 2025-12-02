{
  description = "Fortran and C++ interoperability project";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };
  outputs = { self, nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
          };
        };
      in
      {
        devShell = pkgs.mkShell {
          buildInputs = with pkgs; [
            gcc
            gfortran
            cmake
          ];
          shellHook = ''
            echo "Welcome to the Fortran and C++ interoperability project!"
            echo "To compile, use 'make'"
            echo "To run, use './main'"

            export NIX_BUILD=1
          '';
        };
      }
    );
}

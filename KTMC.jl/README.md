# KTMC
Quailifying the difference between the Omnes and KT

## Content
 - `omega.jl`: example of the unbinned likelihood test
 - `Project.toml`: explicit dependences
 - `Manifest.toml`: all dependences

## To run the notebooks:
1. Clone the repository
```julia
git clone https://github.com/mmikhasenko/KTMC
```
2. get Julia binaries (1.7 recommended)
https://julialang.org/downloads/
3. start julia terminal
4. navigate to the project folder
```julia
cd("Documents/KTMC") # use your path
```
5. Activalte project
```julia
] activate . # `]` to enter the package-manager mode, `backspace` to exit
```
6. install all dependences
```julia
] instantiate # the `]` is not needed if already in the package mode
```
7. Start Pluto server
```julia
using Pluto
Pluto.run()
```
8. Open the notebook, it runs automatically

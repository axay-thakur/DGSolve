# This script automatically installs a list of specified Julia packages
# into the environment of the directory from which it is run.

# --- 1. Define the List of Packages to Install ---
const PACKAGES_TO_INSTALL = [
    "LinearAlgebra",
    "Reexport",
    "Plots",
    "Polynomials",
    "SpecialFunctions",
    "Statistics",
    "Revise"
]

# --- 2. Define the Project Directory ---
const PROJECT_DIR = @__DIR__

println("--- Starting automatic package installation for project in: $(PROJECT_DIR) ---")

try
    using Pkg # Ensure Pkg module is loaded
    
    # Activate the environment. If Project.toml doesn't exist, it will be created.
    Pkg.activate(PROJECT_DIR)
    println("Project environment activated: $(Pkg.project().name).")

    # Add (and install) the specified packages.
    # Pkg.add() implicitly handles resolution, downloading, and updating Project.toml/Manifest.toml.
    println("\nInstalling/updating packages: $(join(PACKAGES_TO_INSTALL, ", "))...")
    Pkg.add(PACKAGES_TO_INSTALL)
    println("All specified packages installed successfully!")
    
    println("\n--- Setup complete! ---")
    println("To use, navigate to $(PROJECT_DIR) in your terminal and run:")
    println("  julia --project=.")
    println("Then, in the Julia REPL, use 'using YourPackage' (e.g., 'using Plots').")

catch e
    println(stderr, "An error occurred during package setup: ", e)
    println(stderr, "Please ensure you have an active internet connection and correct package names.")
    rethrow(e) # Re-throw the error to stop script execution
end
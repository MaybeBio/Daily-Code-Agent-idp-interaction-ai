# General guidelines
- Always make your modifications as minimal and as readable as possible while fulfilling the required objectives in the instructions. Try to avoid any unnecessary code bloat or complexity wherever you can 
- Aim to address a small upstream fix rather than creating an entire separate function with duplicated logic just to address a minor problem
- Try to understand and follow the same design patterns, choices, workflows as in the existing code. Don't depart too far from the types of design decisions and logic that exist already in the code, unless absolutely necessary.

# Modifying the hamiltonian to work with polymonomers
Now it's time to modify the hamiltonian.jl file to implement the energy calculations required for moving and interaction, potentially between all types of monomers, whether Poly or Bulk.

- In principle, the main thing to change is the way the ham and deltaHam functions access the intPot and chemPot values. Change all functions so that they read these values off from the underlying PolyMonomers.

- Once this is done, the intPot and chemPot fields are no longer necessary in Polymer. Remove them from the type and everywhere else they are used, making the necessary modifications to the code.
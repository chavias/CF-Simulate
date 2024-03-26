# CF Simulate

CF Simulate is a curated collection of examples showcasing the application of Continuous Floquet theory up to second order for effective Hamiltonian simulations, alongside exact simulations for comparative analysis, as described in the publication by Chavez et al.

## Features

- **Effective Hamiltonian Simulation**: Simulation of effective Hamiltonians up to second order for various systems.
- **Numerical Simulation**: Perform numerical simulations for various systems.
- **Versatile Toolkit**: Provides a range of tools and functions for similar exact and effective calculations.

## Installation

To install CF Simulate, simply clone the repository:

```bash
git clone https://github.com/yourusername/CF-simulate.git
```

*Note*: This code is designed to be used with SLURM Workload Manager.

## Usage

1. Effective Hamiltonian Simulation for a specific system:

To run the effective Hamiltonian simulation first navigate to the corresponding directory
(e.g C7/effective-simulation/) and edit the Run_effective.sh script to set the simulation parameter and output directory.

```bash
# Navigate to the corresponding directory
cd path/to/directory

# Edit the Run_effective.sh script
nano Run_effective.sh

# Execute the script
./Run_effective.sh
```


2. Numerical Simulation for a specific system:

```bash
# Navigate to the corresponding directory
cd path/to/directory

# Edit the run_sim_time.py script
nano run_sim_time.py

# Edit the gamma_time.csh script
nano gamma_time.csh

# Execute the script
./Run_exact.sh
```

Execute this commands to perform numerical simulations up to second-order accuracy.

<!--Documentation
For detailed instructions and documentation on how to use CF Simulate, please refer to the Documentation file.-->

<!--
## Contribution
Contributions are welcome! If you'd like to contribute to CF Simulate, please follow these steps:
    Fork the repository
    Create your feature branch (git checkout -b feature/YourFeature)
    Commit your changes (git commit -am 'Add some feature')
    Push to the branch (git push origin feature/YourFeature)
    Create a new Pull Request
-->

## License

This project is licensed under the MIT License.

<!-- ## Contact
For any inquiries or suggestions, please feel free to reach out to Your Name.

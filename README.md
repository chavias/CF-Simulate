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

### Effective Hamiltonian Simulation:

- Navigate to the corresponding directory:
      
  ```bash
  cd path/to/directory
  ```
  
- Edit the Run_effective.sh script

  ```bash
  nano Run_effective.sh 
  ```
      
  In this script, you can specify the simulation parameters and output directory as needed.

- Execute the Run_effective.sh script
      
  ```bash
  ./Run_effective.sh
  ```


### Exact Simulation:

- Navigate to the corresponding directory:
      
  ```bash
  cd path/to/directory
  ```

- Compile OPTIM_externshap_time.cc using [GAMMA](https://github.com/tesch1/GAMMA)

  ```bash
  gamma OPTIM_externshap_time.cc -o OPTIM_externshap_time
  ```
      
  This script simulates multispin systems under arbitrary rf irradiation and magic angle spinning.
  The specification of the spin system has to be provided with a .sys file.
  The rf irradiation is specified in a .dat file and set within the run_sim_time.py.
      
- Edit the run_sim_time.py script

  ```bash
  nano run_sim_time.py 
  ```
      
  In this script, you can specify the various simulation parameters.

- Edit the gamma_time.csh script

  ```bash
  nano gamma_time.csh 
  ```
      
  In this script, you can specify the output directory and simulation parameters.

- Execute the Run_effective.sh script
      
  ```bash
  ./Run_effective.sh
  ```

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

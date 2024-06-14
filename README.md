# Multiplexed HGP Codes on Erasure Channel
This is a Python implementation of multiplexed HGP codes simulator.
  
## Simulation Flow
You can simulate the multiplexed quantum communication with HGP code. This implementation use part of [Pruned-Peeling-and-VH-Decoder](https://github.com/Nicholas-Connolly/Pruned-Peeling-and-VH-Decoder) repository included with ["Fast erasure decoder for a class of quantum LDPC codes"](https://arxiv.org/abs/2208.01002) with the permission of the contributors.

Simulation flows are in [this jupyter notebook](https://github.com/parton-quark/Multiplexed_HGP/blob/main/simulation_flow.ipynb).

## folders
- `input_matrices` includes the classical parity check matrices used to generate HGP codes for our numerical simulations.
- `results` includes the raw data for the simulation. 
- `scripts` includes Slurm scripts for running multiple shots of simulations

## Requirements (modules)
- bidict
- json
- numpy
- itertools 

## Citation 
ArXiv submission is available at [HERE](https://arxiv.org/abs/2406.08832). For the citation of this work, please use this bibtex file.

```
@misc{nishio2024multiplexed,
      title={Multiplexed Quantum Communication with Surface and Hypergraph Product Codes}, 
      author={Shin Nishio and Nicholas Connolly and Nicolò Lo Piparo and William John Munro and Thomas Rowan Scruby and Kae Nemoto},
      year={2024},
      eprint={2406.08832},
      archivePrefix={arXiv},
      primaryClass={id='quant-ph' full_name='Quantum Physics' is_active=True alt_name=None in_archive='quant-ph' is_general=False description=None}
}
```

If you have any question on this implementation, send email to parton(at)nii.ac.jp (Shin Nishio).
# FPGA Implementation of Wideband Spectrum Sensing using OMP

Implementation of an Orthogonal Matching Pursuit (OMP) based Wideband Spectrum Sensing (WSS) system in Verilog HDL for FPGA platforms. The design detects sparse occupied frequency bands from compressed measurements using compressed sensing techniques.

---

## Features

- Fully synthesizable Verilog HDL implementation
- OMP-based sparse signal recovery
- Fixed-point arithmetic architecture
- FPGA resource optimized design
- Vivado/XSim simulation support
- MATLAB reference verification

---

## System Architecture

The hardware design consists of five major modules:

1. Measurement Matrix Multiplication  
2. Frobenius Norm Computation  
3. Correlation Energy Estimation  
4. Argmax Detection  
5. Orthogonal Residual Update  

---

## Technologies Used

- Verilog HDL
- Xilinx Vivado
- XSim Simulator
- MATLAB

---

## Optimization Techniques

- Shared MAC-based sequential computation
- Fixed-point integer arithmetic
- BRAM-oriented architecture
- Sequential argmax implementation
- Modified Gram-Schmidt (MGS) orthogonalization

---

## Simulation Results

- Accurate dominant band detection
- Verified against MATLAB reference output
- FPGA synthesis compatible
- Resource-efficient optimized architecture

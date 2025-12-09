# Implementation of Three-Halves Garbled Circuits

This repository contains a high-performance C++ implementation of the **Three-Halves** garbled circuit scheme, built upon the **EMP-toolkit** framework. It features a **branchless** implementation for improved efficiency.


## 📥 Installation

You can download the repository using the following command:

```bash
git clone [https://github.com/GarbleCoder/Implementation-of-threehalves.git](https://github.com/GarbleCoder/Implementation-of-threehalves.git)
cd Implementation-of-threehalves



## 📦 Compilation & Usage

Please execute the following commands in the project root directory to compile and run the test:

1.  **Create and enter the build directory**
    ```bash
    mkdir build && cd build
    ```

2.  **Configure CMake**
    ```bash
    cmake ..
    ```

3.  **Build**
    ```bash
    make
    ```

4.  **Return to the root directory**
    ```bash
    cd ..
    ```

5.  **Run the test**
    ```bash
    ./run ./build/bin/test_garble
    ```

*(Note: Please ensure the `run` script has executable permissions. If not, run `chmod +x run`)*
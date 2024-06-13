# CSE305-FFT

Final project for the course CSE305 - Concurrent and Distributed Computing. In this project, we present implementations and applications of the *Discrete Fourier Transform (DFT)*, focusing on exploring parallel implementation of the approaches. Our methods are based on the *Cooley-Tukey Fast Fourier Transform (FFT)* algorithm, a sequential algorithm that provides $O(n \log{n})$ computation time complexity, which we expand to provide efficient parallel algorithms for the DFT. 

We also implement two important applications of the DFT: Data Compression and Polynomial Multiplication. These serve to show the efficacy of our algorithm, but also to compare the performance of Fourier transform-based approaches against traditional methods. Additionally, we implement the *Number-theoretic Transform (NTT)* with analogous FFT algorithms, which allows us to accurately compute integer polynomial multiplication.

The datasets we used are the [Max Planck Weather Dataset](https://www.kaggle.com/datasets/arashnic/max-planck-weather-dataset) and the [Daily Climate time series dataset](https://www.kaggle.com/datasets/sumanthvrao/daily-climate-time-series-data).

We provide several testing functions in the file `main.cpp`.

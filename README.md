FFT: Fast Fourier Transform

La transformada rápida de fourier es un algoritmo que permite calcular la transformada discreta de fourier propia del espacio L2.

Para calcular la FFT de un conjunto de datos unidimensional use la función fft(signal, n, transform), donde
	- signal es un puntero referente al array al cual se quiere calcular la tranformada
	- n es la longitud del array
	- transform es el puntero que apunta al array donde se almacenara la transformada
Para ver el resultado puede usar la función showTransform(transform,n)

Ejemplo1:
	double signal[4] = {2.0, -9.6, -3.5, 0.6};
	std::complex<double> transform[4];
	fft(signal, 4, transform);
	showTransform(transform, 4);

Ejemplo2:
	std::complex<double> signal[4] = {
		{3.6, 5.6}, {-6.3, 0.5}, {7.0, 9.3}, {6.4, -8.5}
	};
	std::complex<double> transform[4];
	fft(signal, 4, transform);
	showTransform(transform, 4);

También se puede hacer lo mismo para datos float

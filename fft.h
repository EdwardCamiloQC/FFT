#ifndef _FFT_H_
    #define _FFT_H_

    #include <complex>
    #include <iostream>

    void calculateFFT(const double* signal, unsigned int length, std::complex<double>* transform);
    void calculateFFT(const std::complex<double>* signal, unsigned int length, std::complex<double>* transform);
    void calculateFFT(const float* signal, unsigned int length, std::complex<float>* transform);
    void calculateFFT(const std::complex<float>* signal, unsigned int length, std::complex<float>* transform);
    void seeTransform(const std::complex<double>* transform, unsigned int n);
    void seeTransform(const std::complex<float>* transform, unsigned int n);
    void calculateModule(const std::complex<double>* transform, double* module);

    template<typename T1, typename T2>
    void fft(const T1* signal, unsigned int n, T2* transform){
        double log2n = log(n)/log(2);
        double nInt = log2n - floor(log2n);
        if(nInt==0.0){
            calculeFFT(signal, n, transform);
        }else{
            std::cout << "La longitud debe ser 2^n" << std::endl;
        }
    }

    template<typename T1, typename T2>
    void oddValues(const T1* values, unsigned int len, T2* odd){
        if(!(len%2)){
            for(size_t i=0; i<(len/2); i++){
                *odd = *values;
                if(i!=(len/2-1)){
                    values+=2;
                    odd++;
                }else{
                    values-=(i*2);
                    odd-=i;
                }
            }
        }
    }

    template<typename T1, typename T2>
    void evenValues(const T1* values, unsigned int len, T2* even){
        if(!(len%2)){
            values++;
            for(size_t i=0; i<(len/2); i++){
                *even = *values;
                if(i!=(len/2-1)){
                    values+=2;
                    even++;
                }else{
                    values-=(i*2+1);
                    even-=i;
                }
            }
        }
    }
    
#endif

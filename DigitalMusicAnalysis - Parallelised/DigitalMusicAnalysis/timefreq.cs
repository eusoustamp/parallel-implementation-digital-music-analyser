using System;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Collections.Generic;

//Pattarachai Roongsritong
//n10548467
//Parallelisation Project CAB401 QUT

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        private Complex[] twiddles;
        private Complex[] newX;
        private float[][] Y;
        private int N;
        private float fftMax;
        private int array_length;
        private int chunk_size;
        private Stopwatch stopWatch = new Stopwatch();
        public static double stft_time;

        public timefreq(float[] x, int windowSamp)
        {
            int ii;
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];
            Parallel.For(0, wSamp, MainWindow.parallel_options, ii =>
                 {
                     double a = 2 * pi * ii / (double)wSamp;
                     twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
                 });


            timeFreqData = new float[wSamp / 2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];
            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }


            int cols = 2 * nearest / wSamp;

            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }
            stopWatch.Start();
            timeFreqData = stft(compX, wSamp);
            stopWatch.Stop();
            stft_time = stopWatch.Elapsed.TotalMilliseconds;
        }

        float[][] stft(Complex[] x, int wSamp)
        {
            N = x.Length;
            fftMax = 0;
            newX = x;
            array_length = 2 * (int)Math.Floor(N / (double)wSamp)-1;
            chunk_size = (array_length + MainWindow.NUM_THREADS - 1) / MainWindow.NUM_THREADS;

            Y = new float[wSamp / 2][];

            

            Parallel.For(0, wSamp / 2, new ParallelOptions { MaxDegreeOfParallelism = MainWindow.NUM_THREADS }, ll =>
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            });

            

            Thread[] fft_thread = new Thread[MainWindow.NUM_THREADS];

            for (int i = 0; i < MainWindow.NUM_THREADS; i++)
            {
                fft_thread[i] = new Thread(stftOnThread);
                fft_thread[i].Start(i);
            }
            for (int j = 0; j < MainWindow.NUM_THREADS; j++)
            {
                fft_thread[j].Join();
            }

           

            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }
         

            return Y;
        }

        public void stftOnThread(object threadID)
        {
            int id = (int)threadID;
            int start = id * chunk_size;
           
            Complex[] temp = new Complex[wSamp];
            Complex[] tempFFT = new Complex[wSamp];

            for (int ii = start; ii < Math.Min(start + chunk_size, array_length - 1); ii++)
            {
                for (int jj = 0; jj < wSamp; jj++)
                {
                    temp[jj] = newX[ii * (wSamp / 2) + jj];
                }

                tempFFT = fft(temp);

                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

                    if (Y[kk][ii] > fftMax)
                    {
                        fftMax = Y[kk][ii];
                    }
                }
            }

            Complex[] fft(Complex[] x)
            {
                int ii = 0;
                int kk = 0;
                int N = x.Length;

                Complex[] Y = new Complex[N];

                // NEED TO MEMSET TO ZERO?

                if (N == 1)
                {
                    Y[0] = x[0];
                }
                else
                {

                    Complex[] E = new Complex[N / 2];
                    Complex[] O = new Complex[N / 2];
                    Complex[] even = new Complex[N / 2];
                    Complex[] odd = new Complex[N / 2];

                    for (ii = 0; ii < N; ii++)
                    {

                        if (ii % 2 == 0)
                        {
                            even[ii / 2] = x[ii];
                        }
                        if (ii % 2 == 1)
                        {
                            odd[(ii - 1) / 2] = x[ii];
                        }
                    }

                    E = fft(even);
                    O = fft(odd);

                    for (kk = 0; kk < N; kk++)
                    {
                        Y[kk] = E[(kk % (N / 2))] + O[(kk % (N / 2))] * twiddles[kk * wSamp / N];
                    }
                }

                return Y;
            }

        }
    }
}

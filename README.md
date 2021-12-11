# phase-averaged-waveform
MATLAB code for the phase-analysis of continuous data, including EEG and background unit activity (BUA).



`K_PhaseWave`, `K_plotCircPhaseHist_one` was used to create the **phase-averaged waveform** plots in **Figures 7–9** of [Nakamura et al., (2021)](https://doi.org/10.1523/JNEUROSCI.1753-21.2021) for the phase analysis of continuous BUA data. 





# Requirements

The following additional MATLAB programs are required.

- [**CircStat for MATLAB**](https://github.com/circstat/circstat-matlab) by Philipp Berens
    - Specifically the following four functions are required
        - `circ_confmean.m`
        - `circ_mean.m`
        - `circ_r.m`
        - `circ_std.m`
- [**Red Blue Colormap**](https://uk.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap) by Adam Auton
    - `redblue.m`

# K_PhaseWave_demo.mlx

This MATLAB Live Script









## K_PhaseWave.m



```
K_PhaseWave is similar to K_PhaseHist, but works for a pair of ECoG and
LFP waveform signals. It returns values for plotting "phase-triggered
average waveforms".
 
[results, handles] = K_PhaseWave(lfpwaveform,eegwaveform,sourceRate,newRate,b,a,varargin)

INPUT ARGUMENTS
    lfpwaveform  vector of source data from LFP channel

    eegwaveform  vector of source data from EEG channel (lfpwaveform &
                 eegwaveform must have same length)

    sourceRate   sampling rate [Hz] of the input data event and eeg

    newRate      the new sampling rate [Hz] after resample
               In many cases, 1024 is good.

    b, a         b and a as coefficients of a filter transfer function
               b, a must be calculated for newRate rather than souceRate
               b for numeraotr, a for denominator. You can get b and a by:

               [b, a] = butter(n, Wn)

               where n is the order of the Butterworth filter (lowpass), or
               half the order(bandpass). Wn is normalized cuttoff
               frequency, i.e. [cycles per sec] devided by Niquist
               frequency [newRate/2].

               Wn = frequencyHerz/(samplingrateHerz/2)

               The following command can check the stablity of the filter
               fvtool(b,a,'FrequencyScale','log', 'Analysis','info');


OPTIONAL PARAMETER/VALUE PAIRS

    'PlotLinear'   true | false (default)

    'PlotCirc'     true | false (default)

    'Histbin'      number of histogram bins (default = 72)

    'HistType'     'line' (default) | 'bar'

    'Randomization' 
                 'bootstrap' (default) | 'circshift' | 'none'


OUTPUT ARGUMENTS
    results       Structure conttaining following fields

      binmean       Mean per bin
      binstd        SD per bin
      binsem        SEM per bin
      axrad         Phase axis (-pi to pi)

      meanvec       Non-scalar structure with the following fields

          vec       Vector represented as a complex number

                      vec = mean(xy)*2, 

                    where, xy = rad2cmp(x).*y, x is instantaneous phase (in
                    radian), and y is the value of the waveform data
                    (typically in mV), Y = rad2cmp(X) is a function defined
                    by:

                      y = exp(1i * x);

                    where 1i is equivalent to 1*i, and i is the imaginary
                    unit.


          length   abs(vec)
          radian   angle(vec)
          degree   rad2deg(angle(vec))

          bootstrap, circshift    
                   These are results of shuffling by different methods
          p_lessthan 
                   Estimated range of P value
          length   
          radian
          degree

    handles        Structure of graphic handles



NOTE
What about bias correction with ecdf? How does it apply to LFP ddata?

Because here we do not really create histograms of events for Rayleigh
test, but rather create average waveforms in relation to instantaneous
phase in essense, even if the phases are non-uniform, average will be
computed for each bin without bias.
```









# References

- Nakamura KC, Sharott A, Tanaka T, Magill PJ (2021) Input zone-selective dysrhythmia in motor thalamus after dopamine depletion. J Neurosci, ***in press***, https://doi.org/10.1523/JNEUROSCI.1753-21.2021
- Berens P (2009) CircStat: A MATLAB toolbox for circular statistics. J Stat Softw 31:1–21, http://www.jstatsoft.org/v31/i10



# Contacts

Dr Kouichi C. Nakamura

MRC Brain Network Dynamics Unit, University of Oxford

kouichi.c.nakamura@gmail.com

kouichi.nakamura@ndcn.ox.ac.uk




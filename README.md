# Athena_signal_aec_matlab
Matlab version of Athena Signal AEC(subband echo cancellation)

## Aims
- understand WOLA implementation of subband compose/decompose;
- understand subband adaptive filter;

## Simplify
- nlms instead of ipnlms;
- some parameters like filter num, smooth factor etc;
- remove nlp, dtd status etc, keep forefilter/backfilter and copy logic. 

## Reference
- [athena-signal](https://github.com/athena-team/athena-signal)
-  Multirate digital signal processing,  Crochiere, R. E. , & Rabiner, L. R., Chapter 7

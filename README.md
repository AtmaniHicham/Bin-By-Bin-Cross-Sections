```
How to run:
```
-  ./make or python CrossSecBB.py  "Wplusenu"   "plusenu"   "5TeV"  "W^{+} \rightarrow e^{+} \nu, 5.02 TeV" "el" "pT"
  
```
What it does:
```
This framework is used for:
-  Cross sections calculation
-  uncertainties propogation using reconstructed distributions

\begin{equation}
\sigma_\mathrm{fid}=\frac{N^\mathrm{data}-N^\mathrm{bg}}{\mathcal{L} \cdot C_{v}}
\end{equation}
where
\begin{itemize}
    \item for a given channel, $N^\mathrm{data}$ and $N^\mathrm{bg}$ represent the
          number of events of data in the phase space defined in the section, and
          the expected number of background.
    \item $C_{v}$ is a correction factor calculated using simulation, corresponding to
          the ratio of the number of selected events at the detector level and
          the number of events at the particle level in the fiducial phase-space.
          This correction factor allows to correct the observed difference between
          data and simulation (reconstruction, identification, isolation, and trigger).
    \item $\mathcal{L}$ is the integrated luminosity of data.
\end{itemize}
```
Output examples:
```

A three transmon array, where the middle transmon is specifically detuned from the edge transmons, can be used as a tunable coupler between two transmons (edge ones). The edge transmons have native capacitive coupling. From the middle transmon, the edge transmons feel effective coupling with strength $g_1g_2 / \Delta$, where $g_1$ and $g_2$ are the couplings between the edge transmons and the middle transmon, and $\Delta = w - w_c$ is the detuning between the middle transmon and edge transmons . Here I assume that the edge transmons have the same frequency, which is not necessary, and is not the case in the numerical simulation shown below. Then the coupling between the edge transmons is of form $g_1g_2 / \Delta + g$, so with choosing $w_c > w$ we get a negative contribution to the coupling from the middle transmon, and then specifically choosing the frequency we can cancel the coupling of the transmons.

When the coupler is in tune, it cancels the interaction between the edge transmons, both for single and two excitations. Detuning the middle transmon out, we have a effectively two transmon array, which are in resonance. This would not work as a LRU, since the qubit is not protected. The point here is that we need to have disorder on the edge tranmons.

With disorder, we can protect the qubit as we saw in our LRU. When the coupler is in tune, we still get no evolution. But since when $w_c < w$ the contribution is positive, we can have an intermediate value of the coupler frequency where the second excitation coupling is boosted, while we still have good decoupling of the qubit due to the disorder.

![[CouplerInTune.png]]
When the coupler is tuned to off, without also tuning the anharmonicity the coupler will be detuned both in single and two excitations, so that both interactions are turned off.

![[NoCoupler.png]]
The the coupler is completely detuned out (0 Hz, inf Hz), we have effectively a system only consisting of edge transmons. The qubit is protected by disorder on the edge sites.

![[SlightTuning.png]]
With slight tuning in the coupler, we can enhance leakage mobility. This does not increase qubit interaction too much.

![[SlightTuning2.png]]![[SlightTuningC.png]]
Some of the population goes in to the coupler.

The tunable coupler is not a problem for our LRU. At worse it could be just detuned out, removed. We can use it to slightly improve leakage mobility. The qubit is well protected by disorder.

* *Comparable disorder is used in the IQM processor. Why does it work, if here we have too much disorder for the qubits to couple?*
#### Leakage Level In Resonance
In the previous plot, the frequency of the second excited level $|2\rangle$ is not matched to be in resonance. In the figures below it is. Surprisingly, we are able to turn off qubit interaction with the tunable coupler, while having leakage mobility. The caveat is that for realistic anharmonicity, we need to have have relatively low tunable coupler coupling $J_c$, due to which the single excitation spends a lot of time in the coupler. If the coupler is in full resonance in the second level, we get rapid oscillation, so it might be beneficial to have slight detuning for the anharmonicity of the coupler. Note that here the edge transmons have same transition frequency (they are in resonance as qubits and for leakage). 
###### High Coupling $J_c/2\pi = 50$MHz -> $U_c/2\pi = 1140$ MHz   
![[Data/Coupler/FLevelResonance/high_coupling/combined.png]]![[Data/Coupler/FLevelResonance/high_coupling/n2_exact.png]]
###### Low Coupling $J_c/2\pi = 15.8$MHz -> $U_c/2\pi = 240$ MHz   
![[Data/Coupler/FLevelResonance/low_coupling/combined.png]]![[Data/Coupler/FLevelResonance/low_coupling/n2_exact.png]]
With the second level in resonance, it is possible to have leakage mobility while turning off the qubit interaction. But, the coupling strengths dictate the amount of detuning needed to turn of the interaction. This decides the anharmonicity of the coupler. For high coupling strength we would need anharmonicity difference of $\Delta U \approx 1 GHz$, which is unrealistic. We can lower the coupling so that we are not forced to have such large anharmonicity. But for lower coupling a lot of the qubit population will go to the coupler.
* *Is there some intermediate value where we have realistic coupling but tolerable qubit interaction with the coupler?*

Another challenge is that the anharmonicity of the coupler needs to be very specifically engineered. The situation would be that first the edge transmons are engineered, their coupling is decided, and then the coupler needs to be manufactured so that its coupling strength is specific, and anharmonicity is specific.
* *Do superconducting elements exist, which have tunable anharmonicity and are compatible with transmons?*
 
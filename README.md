# MCast: High-Quality Linear Video Transmission with Time and Frequency Diversities
# ParCast: Soft Video Delivery in MIMO-OFDM WLANs

### The Matlab code is for two different joint source and channel coding based wireless transmission systems (SoftCast & ParCast) demonstrated in the result secion in the paper titled: MCast: High-Quality Linear Video Transmission with Time and Frequency Diversities.



In the Parcast+joint and Parcast+average systems, the channel assignment and power allocation at each time slot are the same as those employed by Parcast. Speciﬁcally, at each time slot, both the channels and blocks are sorted in a descend order according to the channel powers and block energy, respectively. Then, each channel is assigned to the corresponding block with the same sorted order, i.e., higher gain channel is assigned to more important block. 

In Softcast+joint and Softcast+average systems, there is no channel assignment procedure and the power allocation at each time slot is the same as that in Softcast. 

Where “average” means that the receiver decodes the received signal of each time slot independently and average the results while “joint” means that the receiver jointly decodes the received signals of all time slots.

## Main features:

### -This contribution aims to compare the performance of different systems, ParCast & SoftCast, as shown in Figure(6) in "MCast: High-Quality Linear Video Transmission with Time and Frequency Diversities".
### -The Rayleigh fading channel model is assumed in this simulation in addition to the presence of Additive White Gaussian Noise (AWGN), assuming the Signal to Noise Ratio (SNR) = 5 dB.
### -Peak Signal to Noise Ratio (PSNR) is used as the standard metric of video quality.
### -The transmission is performed across different time slots = [1,2,3,4,5].



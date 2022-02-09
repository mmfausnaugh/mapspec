import numpy as np
import matplotlib.pyplot as plt

from mapspec.mapspec.spectrum import TextSpec


speclist = ['mcg0811_001_140110.txt',
            'mcg0811_002_140110.txt',
            'mcg0811_003_140110.txt',
            'mcg0811_001_140417.txt',
            'mcg0811_002_140417.txt',
            'mcg0811_003_140417.txt',]

F,(ax1,ax2,ax3) = plt.subplots(3,1,sharex='col')

for spec in speclist:
    s = TextSpec(spec)
    ax1.plot(s.wv, s.f,'-',label=spec)

    s2 = TextSpec('scale_'+spec)
    ax2.plot(s2.wv, s2.f,'-',label='scale_'+spec)

    s3 = TextSpec('scale.h._'+spec)
    ax3.plot(s3.wv, s3.f,'-',label='scale.h._'+spec)


    
F.set_size_inches(12,10)
F.tight_layout()

ax3.set_xlabel('Wavelength ($\AA$)')
ax2.set_ylabel('Flux')

ax1.set_title('Raw Spectra')
ax2.set_title('Rescaling with Gaussian Smoothing')
ax3.set_title('Rescaling with Gauss-Hermite Smoothing')

ax1.legend(loc='upper left')
ax2.legend(loc='upper left')
ax3.legend(loc='upper left')

plt.show()

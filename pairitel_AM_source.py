# -*- coding: utf-8 -*-

# plot Ann Marie's source:
ra_am = 100.25214
dec_am = 9.4877639
i_am = np.where( (np.abs(mycloud['ra'] - ra_am) < 0.001) & (np.abs(mycloud['dec'] - dec_am) < 0.01) )[0]

for b in ['J', 'H', 'K']:
    stuff = zip(mycloud[117].lclist[0]['t'+b+'cal2'], mycloud[117].lclist[0]['m'+b+'cal2'], mycloud[117].lclist[0]['m'+b+'cal2_error']) # sort by first argument in zip.
    stuff.sort()
    (mycloud[117].lclist[0]['t'+b+'cal2'], mycloud[117].lclist[0]['m'+b+'cal2'], mycloud[117].lclist[0]['m'+b+'cal2_error']) = zip(*stuff)


plt.clf()
plt.errorbar(mycloud[117].lclist[0]['tJcal2'], mycloud[117].lclist[0]['mJcal2'], yerr = mycloud[117].lclist[0]['mJcal2_error'], fmt='.-', capsize=0 )
plt.errorbar(mycloud[117].lclist[0]['tHcal2'], mycloud[117].lclist[0]['mHcal2'], yerr = mycloud[117].lclist[0]['mHcal2_error'], fmt='.-', capsize=0  )
plt.errorbar(mycloud[117].lclist[0]['tKcal2'], mycloud[117].lclist[0]['mKcal2'], yerr = mycloud[117].lclist[0]['mKcal2_error'], fmt='.-', capsize=0  )
plt.axis([55550, 56000, 10.5, 13.5])
#plt.errorbar(mycloud[117].lclist[0]['tKcal2'], mycloud[117].lclist[0]['mKcal2'], yerr = mycloud[117].lclist[0]['mKcal1_error'], fmt='.-' )

ascii.write([mycloud[117].lclist[0]['tJcal2'], mycloud[117].lclist[0]['mJcal2'], mycloud[117].lclist[0]['mJcal2_error']] , names=['tJ', 'magJ', 'magJerror'] )

ascii.write([mycloud[117].lclist[0]['tHcal2'], mycloud[117].lclist[0]['mHcal2'], mycloud[117].lclist[0]['mHcal2_error']] , names=['tH', 'magH', 'magHerror'] )

ascii.write([mycloud[117].lclist[0]['tKcal2'], mycloud[117].lclist[0]['mKcal2'], mycloud[117].lclist[0]['mKcal2_error']] , names=['tK', 'magK', 'magKerror'] )


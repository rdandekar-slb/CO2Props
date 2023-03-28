import sys
import math
import matplotlib.pyplot as plt


def co2_solubility(pressure, temperature):
      ai=[1.163, -16.630, 111.073, -376.589, 524.889]
      bi=[0.965, -0.272,0.0923, -0.1008, 0.0998]
      ci=[1.280,-10.757,52.696,-222.395,462.672]
      a = b = c = 0
      for i in range(5):
            a+=ai[i]*pow(10,-3*i)*pow(temperature,i)
            b+=bi[i]*pow(10,-3*i)*pow(temperature,i)
            c+=ci[i]*pow(10,-3*i)*pow(temperature,i)
      c=c*0.001
      p0=(2/math.pi)*((math.asin(b*b))/(c*(1-((2/math.pi)*math.asin(b*b)))))
      if pressure<p0:
            rsw=a*pressure*(1-b*math.sin(math.pi/2 *((c*pressure)/(c*pressure+1))))
      else:
            rsw0=a*p0*(1-b*b*b)
            m=(math.pi/2)*((c*p0)/(c*p0+1))
            m=math.sin(m)+(m/(c*p0+1)*math.cos(m))
            m=a*(1-b*m)
            rsw=rsw0+m*(pressure-p0)
      return rsw

def co2_density(pressure, temperature):
    if pressure<1100:
        return None
    temperature = (temperature-32.0)/1.8
    a = []
    if pressure < 3000:
        b0 = [-214832.2085348, 475.7146002428, -0.3713900186613,
              0.0001228907393482, -1.466408011784E-08]
        b1 = [11681.16599408, -26.19250287624, 0.02072488876536, -
              0.000006930063746226, 8.338008651366E-10]
        b2 = [-230.2236659392, 0.5215134206837, -0.0004169082831078,
              1.406317206628E-07, -1.704242447194E-11]
        b3 = [1.967428940167, -0.004494511089838,
              0.000003622975674137, -1.230995287169E-09, 1.500878861807E-13]
        b4 = [-0.006184842764145, 0.00001423058795982, -
              1.155050860329E-08, 3.94841742804E-12, -4.838826574173E-16]

    else:
        b0 = [689.7382693936, 0.2213692462613, -0.00005118724890479,
              5.517971126745E-09, -2.184152941323E-13]
        b1 = [2.730479206931, -0.006547268255814,
              0.000002019697017603, -2.415814703211E-10, 1.010703706059E-14]
        b2 = [-0.02254102364542, 0.00005982258882656, -
              2.311332097185E-08, 3.121603486524E-12, -1.406620681883E-16]
        b3 = [-0.004651196146917, 0.000002274997412526, -
              4.079557404679E-10, 3.17127108487E-14, -8.957731136447E-19]
        b4 = [0.00003439702234956, -1.88836133766E-08,
              3.893599641874E-12, -3.560785550401E-16, 1.215810469539E-20]

    
    a = [(b0[i]+b1[i]*temperature+b2[i]*temperature*temperature + b3[i]
          * pow(temperature, 3)+b4[i]*pow(temperature, 4)) for i in range(5)]
#     for x in a:
#          print(x)
    rho = a[0]+a[1]*pressure+a[2]*pressure*pressure + \
        a[3]*pow(pressure, 3)+a[4]*pow(pressure, 4)
    rho *= 0.062428
    return rho





if __name__ == "__main__":
#     pressure = 1100
    temperature = 122
#     pressure=[5850]
    
    co2_rho=[]
    co2_rsw=[]
    pressure=[100+i*100 for i in range(100)]
    for p in pressure:
        co2_rho.append(co2_density(p,temperature))
        co2_rsw.append(co2_solubility(p,temperature))

    for p,rho,rsw in zip(pressure, co2_rho, co2_rsw):
        print(p,rho,rsw)
#     plt.plot(pressure,co2_rsw)
    fig,axs=plt.subplots(2)
    fig.suptitle('CO2 properties')
    axs[0].plot(pressure, co2_rho)
    axs[1].plot(pressure, co2_rsw)
    plt.show()


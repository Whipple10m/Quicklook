      real function gain_fit(x)
      common/pawpar/par(3)
      real x

      gain_fit=par(1)+par(2)*(log(x)-par(3))**2
      return
      end

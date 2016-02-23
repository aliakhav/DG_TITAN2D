#ifndef PROPS
#define PROPS

/* struct holds constants for material properties as well as other constants */

struct MatProps{
  double intfrict; //phi_{int}
  double bedfrict; //phi_{bed}
  double porosity; //v_f
  double mu;       //pore fluid viscosity
  double rho;      //density
  double epsilon;  //scaling value
  double gamma;    //slope limiting stuff
  double LENGTH_SCALE, HEIGHT_SCALE, GRAVITY_SCALE;
  double frict_tiny;
  MatProps(double intfrictin, double bedfrictin,
	   double porosityin, double muin, double rhoin,
	   double epsilonin, double gammain, double frict_tinyin,
	   double lscale, double hscale, double gscale) {
    intfrict = intfrictin;
    bedfrict = bedfrictin;
    porosity = porosityin;
    mu = muin;
    rho = rhoin;
    epsilon = epsilonin;
    gamma = gammain;
    frict_tiny = frict_tinyin;
    LENGTH_SCALE = lscale;
    HEIGHT_SCALE = hscale;
    GRAVITY_SCALE = gscale;
  }

};

#endif

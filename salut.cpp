// Sphère 
Su = vec3(-r*sin(u)*cos(v), r*cos(u)*cos(v), 0);
Sv = vec3(-r*cos(u)*sin(v), -r*sin(u)*sin(v), r*cos(v));
Suu = vec3(-r*cos(u)*cos(v), -r*sin(u)*cos(v), 0);
Svv = vec3(-r*cos(u)*cos(v), -r*sin(u)*cos(v), -r*sin(v));
Suv = vec3(r*sin(u)*sin(v), -r*cos(u)*sin(v), 0);
 
// Paraloïde
Su = vec3(a,0,2*h*u);
Sv = vec3(0,b,-2*h*v);
Suu = vec3(0,0,2*h);
Svv = vec3(0,0,-2*h);
Suv = vec3(0,0,0);

// Caténoïde
Su = vec3(a*sinh(u)*cos(v),a*sinh(u)*sin(v),a);
Sv = vec3(-a*cosh(u)*sin(v),a*cosh(u)*cos(v),0);
Suu = vec3(a*cosh(u)*cos(v),a*cosh(u)*sin(v),0);
Svv = vec3(-a*cosh(u)*cos(v),a*cosh(u)*sin(v),0);
Suv = vec3(-a*sinh(u)*sin(v),a*sinh(u)*cos(v),0);

// Hélicoïde

// Pseudo Sphère
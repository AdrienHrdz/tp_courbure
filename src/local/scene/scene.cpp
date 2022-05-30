

#include <GL/glew.h>

#include "scene.hpp"

#include "../../lib/opengl/glutils.hpp"

#include "../../lib/perlin/perlin.hpp"
#include "../../lib/interface/camera_matrices.hpp"

#include "../interface/myWidgetGL.hpp"

#include <cmath>

#include <string>
#include <sstream>
#include "../../lib/mesh/mesh_io.hpp"

#include "../../lib/3d/mat2.hpp"

using namespace cpe;
vec3 colormap(float x);
vec2 curvature(vec3 Su, vec3 Sv, vec3 Suu, vec3 Suv, vec3 Svv);
float normalisationColormap(float Ks, float min, float max);

void scene::build_surface()
{
    // cylindre();
    //sphere();
    // paraboloide_hyperbolique();
    // catenoide();
    // helicoide();
    pseudo_sphere();

    // calcul courbure
    float const min_u = surface.vertex(0,0).x();
    float const min_v = surface.vertex(0,0).y();
    float const max_u = surface.vertex(surface.size_u()-1,surface.size_v()-1).x();
    float const max_v = surface.vertex(surface.size_u()-1,surface.size_v()-1).y();

    // calcul des pas en u et v
    float const du = (max_u-min_u)/(surface.size_u()-1);
    float const dv = (max_v-min_v)/(surface.size_v()-1);

    std::vector<float> K_s_Array;
    std::vector<float> H_s_Array;
    for(int ku = 1; ku < surface.size_u()-1; ku++)
    {
        for(int kv = 1; kv < surface.size_v()-1; kv++)
        {
            // dérivés premières
            vec3 const Su = (surface.vertex(ku-1,kv)-surface.vertex(ku+1,kv))/(2*du);
            vec3 const Sv = (surface.vertex(ku,kv-1)-surface.vertex(ku,kv+1))/(2*dv);
            // dérivés secondes
            vec3 const Suu = (surface.vertex(ku-1,kv)-2*surface.vertex(ku,kv)+surface.vertex(ku+1,kv))/(pow(du,2));
            vec3 const Svv = (surface.vertex(ku,kv-1)-2*surface.vertex(ku,kv)+surface.vertex(ku,kv+1))/(pow(du,2));
            vec3 const Suv = (surface.vertex(ku+1,kv-1)-surface.vertex(ku+1,kv+1)-surface.vertex(ku-1,kv-1)+surface.vertex(ku-1,kv+1))/(2*du)/(2*dv);

            vec2 courbures = curvature(Su,Sv,Suu,Suv,Svv);

            float const K_s = courbures.x();
            float const H_s = courbures.y();
            K_s_Array.push_back(K_s);
            H_s_Array.push_back(H_s);
        }
    }
    float MAX_gauss_curv = *std::max_element(K_s_Array.begin(), K_s_Array.end());
    float MIN_gauss_curv = *std::min_element(K_s_Array.begin(), K_s_Array.end());
    std::cout<<"MAX = "<<MAX_gauss_curv<<std::endl;
    std::cout<<"min = "<<MIN_gauss_curv<<std::endl;
    //int i = 0;
    for(int ku = 1; ku < surface.size_u()-1; ku++)
    {
        for(int kv = 1; kv < surface.size_v()-1; kv++)
        {
            //std::cout<<"ON NORMALISE"<<std::endl;
            // normalisation color
            float K_s = K_s_Array[(kv+1)*(surface.size_u()-1) + ku];
            float x = normalisationColormap(K_s, MIN_gauss_curv, MAX_gauss_curv);
            std::cout<<"x : "<<x<<std::endl;
            // appliquer la couleur au vertex
            vec3 color = colormap(x);
            surface.color(ku,kv) = color;
        }
    }
}

void scene::build_sphere()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = 0.0f;
    float const u_max = 2.0f*M_PI;
    float const v_min = -0.5f*M_PI;
    float const v_max = 0.5f*M_PI;

    float const r = 0.2f;

    vec3 Su;
    vec3 Sv;
    vec3 Suu;
    vec3 Suv;
    vec3 Svv;
    vec2 curvatures;
    float Ks;
    float Hs;
    vec3 color;
    std::vector<float> gauss_curv;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);

            Su = vec3(-r*sin(u)*cos(v), r*cos(u)*cos(v), 0);
            Sv = vec3(-r*cos(u)*sin(v), -r*sin(u)*sin(v), r*cos(v));
            Suu = vec3(-r*cos(u)*cos(v), -r*sin(u)*cos(v), 0);
            Svv = vec3(-r*cos(u)*cos(v), -r*sin(u)*cos(v), -r*sin(v));
            Suv = vec3(r*sin(u)*sin(v), -r*cos(u)*sin(v), 0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            gauss_curv.push_back(Ks);
        }
    }
    float MAX_gauss_curv = *std::max_element(gauss_curv.begin(), gauss_curv.end());
    float MIN_gauss_curv = *std::min_element(gauss_curv.begin(), gauss_curv.end());
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = r*cos(u)*cos(v);
            float const y = r*sin(u)*cos(v);
            float const z = r*sin(v);

            Su = vec3(-r*sin(u)*cos(v), r*cos(u)*cos(v), 0);
            Sv = vec3(-r*cos(u)*sin(v), -r*sin(u)*sin(v), r*cos(v));
            Suu = vec3(-r*cos(u)*cos(v), -r*sin(u)*cos(v), 0);
            Svv = vec3(-r*cos(u)*cos(v), -r*sin(u)*cos(v), -r*sin(v));
            Suv = vec3(r*sin(u)*sin(v), -r*cos(u)*sin(v), 0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            Ks = normalisationColormap(Ks, MIN_gauss_curv, MAX_gauss_curv);
            color = colormap(Ks);
            surface.vertex(ku,kv) = {x,y,z};
            surface.color(ku,kv) = color;
        }
    }
}

void scene::build_paraboloide_hyperbolique()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = -1.0f;
    float const u_max = 1.0f;
    float const v_min = -1.0f;
    float const v_max = 1.0f;

    float const a = 0.2f;
    float const b = 0.2f;
    float const h = 0.2f;

    vec3 Su;
    vec3 Sv;
    vec3 Suu;
    vec3 Suv;
    vec3 Svv;
    vec2 curvatures;
    float Ks;
    float Hs;
    vec3 color;
    std::vector<float> gauss_curv;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);

            Su = vec3(a,0,2*h*u);
            Sv = vec3(0,b,-2*h*v);
            Suu = vec3(0,0,2*h);
            Svv = vec3(0,0,-2*h);
            Suv = vec3(0,0,0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            gauss_curv.push_back(Ks);
        }
    }
    float MAX_gauss_curv = *std::max_element(gauss_curv.begin(), gauss_curv.end());
    float MIN_gauss_curv = *std::min_element(gauss_curv.begin(), gauss_curv.end());
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = a*u;
            float const y = b*v;
            float const z = h*(u*u - v*v);

            Su = vec3(a,0,2*h*u);
            Sv = vec3(0,b,-2*h*v);
            Suu = vec3(0,0,2*h);
            Svv = vec3(0,0,-2*h);
            Suv = vec3(0,0,0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            Ks = normalisationColormap(Ks, MIN_gauss_curv, MAX_gauss_curv);
            color = colormap(Ks);
            surface.vertex(ku,kv) = {x,y,z};
            surface.color(ku,kv) = color;
        }
    }
}

void scene::build_catenoide()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = -1.0f;
    float const u_max = 1.0f;
    float const v_min = 0.0f;
    float const v_max = 2.0f*M_PI;

    float const a = 0.2f;

    vec3 Su;
    vec3 Sv;
    vec3 Suu;
    vec3 Suv;
    vec3 Svv;
    vec2 curvatures;
    float Ks;
    float Hs;
    vec3 color;
    std::vector<float> gauss_curv;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);

            Su = vec3(a*sinh(u)*cos(v),a*sinh(u)*sin(v),a);
            Sv = vec3(-a*cosh(u)*sin(v),a*cosh(u)*cos(v),0);
            Suu = vec3(a*cosh(u)*cos(v),a*cosh(u)*sin(v),0);
            Svv = vec3(-a*cosh(u)*cos(v),-a*cosh(u)*sin(v),0);
            Suv = vec3(-a*sinh(u)*sin(v),a*sinh(u)*cos(v),0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            gauss_curv.push_back(Ks);
        }
    }
    float MAX_gauss_curv = *std::max_element(gauss_curv.begin(), gauss_curv.end());
    float MIN_gauss_curv = *std::min_element(gauss_curv.begin(), gauss_curv.end());
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = a*cosh(u)*cos(v);
            float const y = a*cosh(u)*sin(v);
            float const z = a*u;

            Su = vec3(a*sinh(u)*cos(v),a*sinh(u)*sin(v),a);
            Sv = vec3(-a*cosh(u)*sin(v),a*cosh(u)*cos(v),0);
            Suu = vec3(a*cosh(u)*cos(v),a*cosh(u)*sin(v),0);
            Svv = vec3(-a*cosh(u)*cos(v),-a*cosh(u)*sin(v),0);
            Suv = vec3(-a*sinh(u)*sin(v),a*sinh(u)*cos(v),0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            Ks = normalisationColormap(Ks, MIN_gauss_curv, MAX_gauss_curv);
            color = colormap(Ks);
            surface.vertex(ku,kv) = {x,y,z};
            surface.color(ku,kv) = color;
        }
    }
}

void scene::build_helicoide_droit()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = 0.0f;
    float const u_max = 2.0f;
    float const v_min = 0.0f;
    float const v_max = 2.0f*M_PI;

    float const a = 0.2f;
    float const h = 0.5f;

    vec3 Su;
    vec3 Sv;
    vec3 Suu;
    vec3 Suv;
    vec3 Svv;
    vec2 curvatures;
    float Ks;
    float Hs;
    vec3 color;
    std::vector<float> gauss_curv;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);

            Su = vec3(a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                      a*cosh(u-v)*sin(u+v) + a*sinh(u-v)*cos(u+v),
                      h);
            Sv = vec3(-a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                      -a*cosh(u-v)*sin(u+v) + a*sinh(u-v)*cos(u+v),
                      h);
            Suu = vec3(a*sinh(u-v)*cos(u+v) - a*cosh(u-v)*sin(u+v) - a*cosh(u-v)*sin(u+v) - a*sinh(u-v)*cos(u+v),
                       a*sinh(u-v)*sin(u+v) + a*cosh(u-v)*cos(u+v) + a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                       0);
            Svv = vec3(a*sinh(u-v)*cos(u+v) + a*cosh(u-v)*sin(u+v) + a*cosh(u-v)*sin(u+v) - a*sinh(u-v)*cos(u+v),
                       a*sinh(u-v)*sin(u+v) - a*cosh(u-v)*cos(u+v) - a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                       0);
            Suv = vec3(-a*sinh(u-v)*cos(u+v) - a*cosh(u-v)*sin(u+v) + a*cosh(u-v)*sin(u+v) - a*sinh(u-v)*cos(u+v),
                       -a*sinh(u-v)*sin(u+v) + a*cosh(u-v)*cos(u+v) - a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                       0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            gauss_curv.push_back(Ks);
        }
    }
    float MAX_gauss_curv = *std::max_element(gauss_curv.begin(), gauss_curv.end());
    float MIN_gauss_curv = *std::min_element(gauss_curv.begin(), gauss_curv.end());
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = u*cos(v);//a*sinh(u-v)*cos(u+v);
            float const y = u*sin(v);//a*sinh(u-v)*sin(u+v);
            float const z = h*v;//h*(u+v);

            Su = vec3(a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                      a*cosh(u-v)*sin(u+v) + a*sinh(u-v)*cos(u+v),
                      h);
            Sv = vec3(-a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                      -a*cosh(u-v)*sin(u+v) + a*sinh(u-v)*cos(u+v),
                      h);
            Suu = vec3(a*sinh(u-v)*cos(u+v) - a*cosh(u-v)*sin(u+v) - a*cosh(u-v)*sin(u+v) - a*sinh(u-v)*cos(u+v),
                       a*sinh(u-v)*sin(u+v) + a*cosh(u-v)*cos(u+v) + a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                       0);
            Svv = vec3(a*sinh(u-v)*cos(u+v) + a*cosh(u-v)*sin(u+v) + a*cosh(u-v)*sin(u+v) - a*sinh(u-v)*cos(u+v),
                       a*sinh(u-v)*sin(u+v) - a*cosh(u-v)*cos(u+v) - a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                       0);
            Suv = vec3(-a*sinh(u-v)*cos(u+v) - a*cosh(u-v)*sin(u+v) + a*cosh(u-v)*sin(u+v) - a*sinh(u-v)*cos(u+v),
                       -a*sinh(u-v)*sin(u+v) + a*cosh(u-v)*cos(u+v) - a*cosh(u-v)*cos(u+v) - a*sinh(u-v)*sin(u+v),
                       0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            Ks = normalisationColormap(Ks, MIN_gauss_curv, MAX_gauss_curv);
            color = colormap(Ks);
            surface.vertex(ku,kv) = {x,y,z};
            surface.color(ku,kv) = color;
        }
    }
}

void scene::build_pseudo_sphere()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = -3.0f;
    float const u_max = 3.0f;
    float const v_min = -1.0f*M_PI;
    float const v_max = 1.0f*M_PI;

    float const a = 0.2f;

    vec3 Su;
    vec3 Sv;
    vec3 Suu;
    vec3 Suv;
    vec3 Svv;
    vec2 curvatures;
    float Ks;
    float Hs;
    vec3 color;
    std::vector<float> gauss_curv;

    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);

            Su = vec3(-a*cos(v)*sinh(u)/(cosh(u)*cosh(u)),
                      -a*sin(v)*sinh(u)/(cosh(u)*cosh(u)),
                      a*tanh(u)*tanh(u));
            Sv = vec3(-a*sin(v)/cosh(u),
                      a*cos(v)/cosh(u),
                      0);
            Suu = vec3(-a*cos(v)*(cosh(u)*cosh(u)*cosh(u) - 2*sinh(u)*sinh(u)*cosh(u))/(cosh(u)*cosh(u)*cosh(u)*cosh(u)),
                       -a*sin(v)*(cosh(u)*cosh(u)*cosh(u) - 2*sinh(u)*sinh(u)*cosh(u))/(cosh(u)*cosh(u)*cosh(u)*cosh(u)),
                       a*(1-tanh(u)*tanh(u))*tanh(u));
            Svv = vec3(-a*cos(v)/cosh(u),
                       -a*sin(v)/cosh(u),
                       0);
            Suv = vec3(a*sinh(u)*sin(v)/(cosh(u)*cosh(u)),
                       -a*sinh(u)*cos(v)/(cosh(u)*cosh(u)),
                       0);

            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            gauss_curv.push_back(Ks);
        }
    }
    float MAX_gauss_curv = *std::max_element(gauss_curv.begin(), gauss_curv.end());
    float MIN_gauss_curv = *std::min_element(gauss_curv.begin(), gauss_curv.end());
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = a*cos(v)/cosh(u);
            float const y = a*sin(v)/cosh(u);
            float const z = a*(u-tanh(u));

//            Su = vec3(-a*cos(v)*sinh(u)/(cosh(u)*cosh(u)),
//                      -a*sin(v)*sinh(u)/(cosh(u)*cosh(u)),
//                      a*tanh(u)*tanh(u));
//            Sv = vec3(-a*sin(v)/cosh(u),
//                      a*cos(v)/cosh(u),
//                      0);
//            Suu = vec3(-a*cos(v)*(cosh(u)*cosh(u)*cosh(u) - 2*sinh(u)*sinh(u)*cosh(u))/(cosh(u)*cosh(u)*cosh(u)*cosh(u)),
//                       -a*sin(v)*(cosh(u)*cosh(u)*cosh(u) - 2*sinh(u)*sinh(u)*cosh(u))/(cosh(u)*cosh(u)*cosh(u)*cosh(u)),
//                       a*(1-tanh(u)*tanh(u))*tanh(u));
//            Svv = vec3(-a*cos(v)/cosh(u),
//                       -a*sin(v)/cosh(u),
//                       0);
//            Suv = vec3(a*sinh(u)*sin(v)/(cosh(u)*cosh(u)),
//                       -a*sinh(u)*cos(v)/(cosh(u)*cosh(u)),
//                       0);
            Su = vec3(-a*tanh(u)*cos(v)/cosh(u),
                      -a*tanh(u)*sin(v)/cosh(u),
                      a*tanh(u)*tanh(u));
            Sv = vec3(-a*sin(v)/cosh(u),
                      -a*cos(v)/cosh(u),
                      0);
            Suu = vec3(0.5*a*(cosh(2*u)-3)*cos(v)/(cosh(u)*cosh(u)*cosh(u)),
                       0.5*a*(cosh(2*u)-3)*sin(v)/(cosh(u)*cosh(u)*cosh(u)),
                       2*a*tanh(u)/(cosh(u)*cosh(u)));
            Svv = vec3(-a*cos(v)/cosh(u),
                       -a*sin(v)/cosh(u),
                       0);
            Suv = vec3(a*tanh(u)*sin(v)/cosh(u),
                       -a*tanh(u)*cos(v)/cosh(u),
                       0);


            curvatures = curvature(Su,Sv,Suu,Suv,Svv);
            Ks = curvatures(0);
            Ks = normalisationColormap(Ks, MIN_gauss_curv, MAX_gauss_curv);
            color = colormap(Ks);
            surface.vertex(ku,kv) = {x,y,z};
            surface.color(ku,kv) = color;
        }
    }
}

void scene::sphere()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = 0.0f;
    float const u_max = 2.0f*M_PI;
    float const v_min = -0.5f*M_PI;
    float const v_max = 0.5f*M_PI;

    float const r = 0.2f;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = r*cos(u)*cos(v);
            float const y = r*sin(u)*cos(v);
            float const z = r*sin(v);

            surface.vertex(ku,kv) = {x,y,z};
        }
    }
}

void scene::paraboloide_hyperbolique()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = -1.0f;
    float const u_max = 1.0f;
    float const v_min = -1.0f;
    float const v_max = 1.0f;

    float const a = 1.0f;
    float const b = 1.0f;
    float const h = 1.0f;

    for(int ku = 0; ku < Nu; ku++)
    {
        for(int kv = 0; kv < Nv; kv++)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);

            float const x = a*u;
            float const y = b*v;
            float const z = h*(pow(u,2)-pow(v,2));

            surface.vertex(ku,kv) = {x,y,z};
        }
    }
}

void scene::catenoide()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = -1.0f;
    float const u_max = 1.0f;
    float const v_min = 0.0f;
    float const v_max = 2.0f*M_PI;

    float const a = 0.2f;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = a*cosh(u)*cos(v);
            float const y = a*cosh(u)*sin(v);
            float const z = a*u;

            surface.vertex(ku,kv) = {x,y,z};
        }
    }
}

void scene::helicoide_droit()
{

}

void scene::pseudo_sphere()
{
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = -3.0f;
    float const u_max = 3.0f;
    float const v_min = -1.0f*M_PI;
    float const v_max = 1.0f*M_PI;

    float const a = 0.2f;


    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = a*cos(v)/cosh(u);
            float const y = a*sin(v)/cosh(u);
            float const z = a*(u-tanh(u));

            surface.vertex(ku,kv) = {x,y,z};
        }
    }
}

void scene::load_scene()
{
    load_common_data();

    // Choix de la surface à tracer
    // build_surface();
    // build_sphere();
    // build_paraboloide_hyperbolique();
    // build_catenoide();
    // build_helicoide_droit();
    build_pseudo_sphere();


    surface.fill_normal();
    surface.fill_empty_field_by_default();

    surface_opengl.fill_vbo(surface);

}

void scene::draw_scene()
{
    //Setup uniform parameters
    glUseProgram(shader_program_id);                                                                           PRINT_OPENGL_ERROR();

    //Get cameras parameters (modelview,projection,normal).
    camera_matrices const& cam=pwidget->camera();

    //Set Uniform data to GPU
    glUniformMatrix4fv(get_uni_loc(shader_program_id,"camera_modelview"),1,false,cam.modelview.pointer());     PRINT_OPENGL_ERROR();
    glUniformMatrix4fv(get_uni_loc(shader_program_id,"camera_projection"),1,false,cam.projection.pointer());   PRINT_OPENGL_ERROR();
    glUniformMatrix4fv(get_uni_loc(shader_program_id,"normal_matrix"),1,false,cam.normal.pointer());           PRINT_OPENGL_ERROR();


    glBindTexture(GL_TEXTURE_2D,texture_default);  PRINT_OPENGL_ERROR();
    surface_opengl.draw();

}


scene::scene()
    :shader_program_id(0)
{}


GLuint scene::load_texture_file(std::string const& filename)
{
    return pwidget->load_texture_file(filename);
}

void scene::set_widget(myWidgetGL* widget_param)
{
    pwidget=widget_param;
}


void scene::load_common_data()
{
    texture_default=load_texture_file("data/white.jpg");
    shader_program_id = read_shader("shaders/shader_mesh.vert",
                                    "shaders/shader_mesh.frag"); PRINT_OPENGL_ERROR();
}

vec3 colormap(float x)
{
    float r = 0.0, g = 0.0, b = 0.0;

    if (x < 0.0) {
        r = 127.0 / 255.0;
    } else if (x <= 1.0 / 9.0) {
        r = 1147.5 * (1.0 / 9.0 - x) / 255.0;
    } else if (x <= 5.0 / 9.0) {
        r = 0.0;
    } else if (x <= 7.0 / 9.0) {
        r = 1147.5 * (x - 5.0 / 9.0) / 255.0;
    } else {
        r = 1.0;
    }

    if (x <= 1.0 / 9.0) {
        g = 0.0;
    } else if (x <= 3.0 / 9.0) {
        g = 1147.5 * (x - 1.0 / 9.0) / 255.0;
    } else if (x <= 7.0 / 9.0) {
        g = 1.0;
    } else if (x <= 1.0) {
        g = 1.0 - 1147.5 * (x - 7.0 / 9.0) / 255.0;
    } else {
        g = 0.0;
    }

    if (x <= 3.0 / 9.0) {
        b = 1.0;
    } else if (x <= 5.0 / 9.0) {
        b = 1.0 - 1147.5 * (x - 3.0 / 9.0) / 255.0;
    } else {
        b = 0.0;
    }

    return vec3(r, g, b);
}

vec2 curvature(vec3 Su, vec3 Sv, vec3 Suu, vec3 Suv, vec3 Svv)
{
    // First fundamental form
    mat2 Is = { dot(Su,Su), dot(Su,Sv),
                dot(Su,Sv), dot(Sv,Sv)};
    // normal
    vec3 n = normalized(cross(Su,Sv));
    // Second fundamental form
    mat2 IIs = {dot(n,Suu), dot(n,Suv),
                dot(n,Suv), dot(n,Svv)};
    // Weingarten matrix
    mat2 Ws = IIs * inverse(Is);
    // eigenvalues are principals curvatures
    vec2 LAMBDA = eigenvalue(Ws);
    float lambda1 = LAMBDA(0);
    float lambda2 = LAMBDA(1);

    // Gauss curvature
    float Ks = lambda1*lambda2;
    // Mean curvature
    float Hs = 0.5*(lambda1+lambda2);
    //std::cout<< "Ks : "<<Ks<<"  ;Hs : "<<Hs<<std::endl;
    return {Ks,Hs};
}

float normalisationColormap(float Ks, float min, float max)
{
//    Pas tout sûr de cette fonction,
//    à voir si on trouve une meilleure façon
//    float len = std::max( abs(min), abs(max));
//    if(Ks<0)
//    {
//        return (len-Ks)/(2*len);
//    }
//    else
//    {
//        return (len+Ks)/(2*len);
//    }
    return (Ks-min)/(max-min);
}

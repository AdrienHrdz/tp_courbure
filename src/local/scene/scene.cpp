

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
    int const Nu=200;
    int const Nv=200;
    surface.set_plane_xy_unit(Nu,Nv);

    float const u_min = 0.0f;
    float const u_max = 2.0f*M_PI;
    float const v_min = 0.0f;
    float const v_max = 2.0f;

    float const r = 0.2f;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = r*cos(u);
            float const y = r*sin(u);
            float const z = v;

            surface.vertex(ku,kv) = {x,y,z};
            surface.color(ku,kv) = {1.0f,v_n,1.0f};
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
    float const u_max = 1.0f;
    float const v_min = 0.0f;
    float const v_max = 4.0f*M_PI;

    float const a = 0.2f;
    float const h = 2.0f;

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
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            float const u_n = static_cast<float>(ku)/(Nu-1);
            float const v_n = static_cast<float>(kv)/(Nv-1);

            float const u = u_min + u_n * (u_max-u_min);
            float const v = v_min + v_n * (v_max-v_min);


            float const x = a*sinh(u-v)*cos(u+v);
            float const y = a*sinh(u-v)*sin(u+v);
            float const z = h*(u+v);

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
            Ks = curvatures(0)/MAX_gauss_curv;
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
            Ks = curvatures(0)/MAX_gauss_curv;
            color = colormap(Ks);
            surface.vertex(ku,kv) = {x,y,z};
            surface.color(ku,kv) = color;
        }
    }
}

void scene::load_scene()
{
    load_common_data();

    // Choix de la surface à tracer
    //build_surface();
    build_sphere();
    //build_paraboloide_hyperbolique();
    //build_catenoide();
    //build_helicoide_droit();
    //build_pseudo_sphere();


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

//vec3 colormap(float x)
/*
{
    if (x < 0.0) {
        return vec3(0.0, 0.0, 0.0);
    } else if (1.0 < x) {
        return vec3(0.0, 0.0, 0.0);
    }
    if (x < 3.1250000000000000e-02) {
        float dx = x - 1.5625000000000000e-02;
        return ((vec3(-1.4151576683620706e+02,  2.4271369358056621e+01,  4.5510373586485706e+01 ) * dx
               + vec3( 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00 )) * dx
               + vec3( 2.6007355728658488e-01,  1.4968553250962457e+00,  3.0913652594248364e+00 )) * dx
               + vec3( 2.0810000000000001e-01,  1.6630000000000000e-01,  5.2920000000000000e-01 );
    } else if (x < 4.6875000000000000e-02) {
        float dx = x - 3.1250000000000000e-02;
        return ((vec3(-5.1390461057291191e+01,  1.2211762733842230e+01, -1.2843448884986955e+01 ) * dx
               + vec3(-6.6335515704472066e+00,  1.1377204386589042e+00,  2.1332987618665173e+00 )) * dx
               + vec3( 1.5642431399834725e-01,  1.5146322069502911e+00,  3.1246980525790007e+00 )) * dx
               + vec3( 2.1162380952380999e-01,  1.8978095238095199e-01,  5.7767619047619101e-01 );
    } else if (x < 6.2500000000000000e-02) {
        float dx = x - 4.6875000000000000e-02;
        return ((vec3(-1.4725107464858192e+02,  1.3014608277362621e+01,  5.8634219534912217e+00 ) * dx
               + vec3(-9.0424794325077311e+00,  1.7101468168077587e+00,  1.5312620953827538e+00 )) * dx
               + vec3(-8.8513670422823654e-02,  1.5591301328169576e+00,  3.1819568159735203e+00 )) * dx
               + vec3( 2.1225238095238100e-01,  2.1377142857142900e-01,  6.2697142857142896e-01 );
    } else if (x < 7.8125000000000000e-02) {
        float dx = x - 6.2500000000000000e-02;
        return ((vec3(-2.1469400225321081e+02, -1.4338005366630648e+01, -4.1817857976177763e+01 ) * dx
               + vec3(-1.5944873556660008e+01,  2.3202065798091316e+00,  1.8061099994526548e+00 )) * dx
               + vec3(-4.7894106087856969e-01,  1.6221044046390967e+00,  3.2341032549553237e+00 )) * dx
               + vec3( 2.0810000000000001e-01,  2.3860000000000001e-01,  6.7708571428571396e-01 );
    } else if (x < 9.3750000000000000e-02) {
        float dx = x - 7.8125000000000000e-02;
        return ((vec3(-2.8846495443400278e+02,  2.0037550842697090e+02,  1.1771734328417965e+02 ) * dx
               + vec3(-2.6008654912279265e+01,  1.6481125782483199e+00, -1.5410209318067788e-01 )) * dx
               + vec3(-1.1344649432057459e+00,  1.6841093914837442e+00,  3.2599158784908235e+00 )) * dx
               + vec3( 1.9590476190476200e-01,  2.6445714285714300e-01,  7.2789999999999999e-01 );
    } else if (x < 1.0937500000000000e-01) {
        float dx = x - 9.3750000000000000e-02;
        return ((vec3(-5.4509738001026233e+02,  5.1696771659011155e+01, -6.5374637230314454e+02 ) * dx
               + vec3(-3.9530449651373146e+01,  1.1040714535762580e+01,  5.3638983732652425e+00 )) * dx
               + vec3(-2.1585134520128149e+00,  1.8823723151401646e+00,  3.3413189453671448e+00 )) * dx
               + vec3( 1.7072857142857101e-01,  2.9193809523809500e-01,  7.7924761904761897e-01 );
    } else if (x < 1.2500000000000000e-01) {
        float dx = x - 1.0937500000000000e-01;
        return ((vec3( 2.3639968744743715e+03, -8.1036503315845437e+02, -8.1573269216733058e+02 ) * dx
               + vec3(-6.5081889339354191e+01,  1.3464000707278728e+01, -2.5280462828444659e+01 )) * dx
               + vec3(-3.7930812487429293e+00,  2.2652584908126849e+00,  3.0301226257549660e+00 )) * dx
               + vec3( 1.2527142857142901e-01,  3.2424285714285700e-01,  8.3027142857142899e-01 );
    } else if (x < 1.4062500000000000e-01) {
        float dx = x - 1.2500000000000000e-01;
        return ((vec3( 1.4125902630655582e+03,  2.5375056097507152e+02,  9.0826266478267496e+02 ) * dx
               + vec3( 4.5730464151631985e+01, -2.4521860222023822e+01, -6.3517932773788282e+01 )) * dx
               + vec3(-4.0954472673010889e+00,  2.0924794358947931e+00,  1.6426476944700765e+00 )) * dx
               + vec3( 5.9133333333333399e-02,  3.5983333333333301e-01,  8.6833333333333296e-01 );
    } else if (x < 1.5625000000000000e-01) {
        float dx = x - 1.4062500000000000e-01;
        return ((vec3(-1.9850459267366693e+03,  1.4738473211499172e+02,  2.4976683303608979e+02 ) * dx
               + vec3( 1.1194563273283002e+02, -1.2627302676317344e+01, -2.0943120362100398e+01 )) * dx
               + vec3(-1.6317582534813697e+00,  1.5120237656082123e+00,  3.2294373922181602e-01 )) * dx
               + vec3( 1.1695238095238101e-02,  3.8750952380952403e-01,  8.8195714285714299e-01 );
    } else if (x < 1.7187500000000000e-01) {
        float dx = x - 1.5625000000000000e-01;
        return ((vec3(-1.3211246088080517e+02,  6.1731462945951478e+01,  9.6199145930320853e+01 ) * dx
               + vec3( 1.8896604917048652e+01, -5.7186433584271068e+00, -9.2353000635336890e+00 )) * dx
               + vec3( 4.1265170979798449e-01,  1.2253683588153301e+00, -1.4859407992871662e-01 )) * dx
               + vec3( 5.9571428571428596e-03,  4.0861428571428599e-01,  8.8284285714285704e-01 );
    } else if (x < 1.8750000000000000e-01) {
        float dx = x - 1.7187500000000000e-01;
        return ((vec3(-2.4276114402580023e+02,  1.8878292291818184e+01,  5.4500811814199913e+01 ) * dx
               + vec3( 1.2703833313260910e+01, -2.8249810328356313e+00, -4.7259650980498993e+00 )) * dx
               + vec3( 9.0640855714657143e-01,  1.0918742277018498e+00, -3.6673884807846019e-01 )) * dx
               + vec3( 1.6514285714285700e-02,  4.2659999999999998e-01,  8.7863333333333304e-01 );
    } else if (x < 2.0312500000000000e-01) {
        float dx = x - 1.8750000000000000e-01;
        return ((vec3(-2.4875702015890445e+02,  2.7531596458333780e+01,  1.1605149669749400e+01 ) * dx
               + vec3( 1.3244046870515243e+00, -1.9400610816566539e+00, -2.1712395442592785e+00 )) * dx
               + vec3( 1.1255997759014531e+00,  1.0174204446629080e+00, -4.7450767061454108e-01 )) * dx
               + vec3( 3.2852380952381001e-02,  4.4304285714285702e-01,  8.7195714285714299e-01 );
    } else if (x < 2.1875000000000000e-01) {
        float dx = x - 2.0312500000000000e-01;
        return ((vec3( 6.6879357994795782e+01,  3.3156266362545779e+00,  3.1398894268734253e+01 ) * dx
               + vec3(-1.0336080632897122e+01, -6.4951749767225808e-01, -1.6272481534897754e+00 )) * dx
               + vec3( 9.8479233924761567e-01,  9.7695827936089374e-01, -5.3385904089187008e-01 )) * dx
               + vec3( 4.9814285714285700e-02,  4.5857142857142902e-01,  8.6405714285714297e-01 );
    } else if (x < 2.3437500000000000e-01) {
        float dx = x - 2.1875000000000000e-01;
        return ((vec3(-3.7807546774099214e+00,  2.9110963663947160e+01,  2.0085673255558202e+01 ) * dx
               + vec3(-7.2011107268910699e+00, -4.9409749909782474e-01, -1.5542498464285720e-01 )) * dx
               + vec3( 7.1077372425092522e-01,  9.5908929503636120e-01, -5.6171330867519242e-01 )) * dx
               + vec3( 6.2933333333333299e-02,  4.7369047619047600e-01,  8.5543809523809500e-01 );
    } else if (x < 2.5000000000000000e-01) {
        float dx = x - 2.3437500000000000e-01;
        return ((vec3(-1.8052110713761824e+01,  7.5676044216235097e+00,  2.6820241280346455e+01 ) * dx
               + vec3(-7.3783336023946600e+00,  8.7047892264969851e-01,  7.8609094921143352e-01 )) * dx
               + vec3( 4.8296990660583561e-01,  9.6497025477935916e-01, -5.5185915297880839e-01 )) * dx
               + vec3( 7.2266666666666701e-02,  4.8866666666666703e-01,  8.4670000000000001e-01 );
    } else if (x < 2.6562500000000000e-01) {
        float dx = x - 2.5000000000000000e-01;
        return ((vec3(-8.5042116753280467e+01,  3.9234694840689350e+01,  6.3623990194130904e+01 ) * dx
               + vec3(-8.2245262921022455e+00,  1.2252103799133005e+00,  2.0432897592276738e+00 )) * dx
               + vec3( 2.3917522075432149e-01,  9.9771540013190607e-01, -5.0765007940944740e-01 )) * dx
               + vec3( 7.7942857142857203e-02,  5.0398571428571404e-01,  8.3837142857142899e-01 );
    } else if (x < 2.8125000000000000e-01) {
        float dx = x - 2.6562500000000000e-01;
        return ((vec3(-4.4981860368289709e+01,  3.5222378119677195e+01,  1.8276940800992332e+01 ) * dx
               + vec3(-1.2210875514912267e+01,  3.0643367005706139e+00,  5.0256642995775600e+00 )) * dx
               + vec3(-8.0127932480280273e-02,  1.0647395732644671e+00, -3.9719767224061564e-01 )) * dx
               + vec3( 7.9347619047619000e-02,  5.2002380952381000e-01,  8.3118095238095202e-01 );
    } else if (x < 2.9687500000000000e-01) {
        float dx = x - 2.8125000000000000e-01;
        return ((vec3( 8.8958586797831074e+01, -6.4031864461777545e+01, -5.4343639113056135e+01 ) * dx
               + vec3(-1.4319400219675847e+01,  4.7153856749304826e+00,  5.8823958996240755e+00 )) * dx
               + vec3(-4.9466349083321959e-01,  1.1862977353816719e+00, -2.2675923162809006e-01 )) * dx
               + vec3( 7.4942857142857103e-02,  5.3754285714285699e-01,  8.2627142857142899e-01 );
    } else if (x < 3.1250000000000000e-01) {
        float dx = x - 2.9687500000000000e-01;
        return ((vec3( 2.3465669412937996e+02, -7.4943148843863256e+01, -1.7040059387215410e+02 ) * dx
               + vec3(-1.0149466463527515e+01,  1.7138920282846606e+00,  3.3350378161995691e+00 )) * dx
               + vec3(-8.7698953275827207e-01,  1.2867551994944084e+00, -8.2736829818345611e-02 )) * dx
               + vec3( 6.4057142857142799e-02,  5.5698571428571397e-01,  8.2395714285714305e-01 );
    } else if (x < 3.2812500000000000e-01) {
        float dx = x - 3.1250000000000000e-01;
        return ((vec3( 3.5054309382746595e+02, -7.5598816353949772e+01, -5.9224118732067950e+01 ) * dx
               + vec3( 8.5006607378717081e-01, -1.7990680737714295e+00, -4.6524900215576546e+00 )) * dx
               + vec3(-1.0222926638479650e+00,  1.2854243237836778e+00, -1.0332202052706571e-01 )) * dx
               + vec3( 4.8771428571428597e-02,  5.7722380952381003e-01,  8.2282857142857202e-01 );
    } else if (x < 3.4375000000000000e-01) {
        float dx = x - 3.2812500000000000e-01;
        return ((vec3(-1.3511844086782639e+02,  2.1571557117596814e+01,  6.5912402293741552e+00 ) * dx
               + vec3( 1.7281773596949638e+01, -5.3427625903628249e+00, -7.4286205871233397e+00 )) * dx
               + vec3(-7.3898266899270237e-01,  1.1738332196565799e+00, -2.9208937378770627e-01 )) * dx
               + vec3( 3.4342857142857203e-02,  5.9658095238095199e-01,  8.1985238095238100e-01 );
    } else if (x < 3.5937500000000000e-01) {
        float dx = x - 3.4375000000000000e-01;
        return ((vec3(-1.6458788273706924e+02,  1.0533768835542057e+01,  3.0362548290707878e+01 ) * dx
               + vec3( 1.0948096681270275e+01, -4.3315958504754741e+00, -7.1196562013714262e+00 )) * dx
               + vec3(-2.9789094589551629e-01,  1.0226713690184817e+00, -5.1940619860793691e-01 )) * dx
               + vec3( 2.6499999999999999e-02,  6.1370000000000002e-01,  8.1350000000000000e-01 );
    } else if (x < 3.7500000000000000e-01) {
        float dx = x - 3.5937500000000000e-01;
        return ((vec3(-1.0406115199344315e+02,  1.9929786587720105e+01,  3.6734795179105028e+01 ) * dx
               + vec3( 3.2330396779701545e+00, -3.8378254363094402e+00, -5.6964117502444944e+00 )) * dx
               + vec3(-7.6310690282384588e-02,  8.9502416141246732e-01, -7.1965726035193567e-01 )) * dx
               + vec3( 2.3890476190476202e-02,  6.2866190476190498e-01,  8.0376190476190501e-01 );
    } else if (x < 3.9062500000000000e-01) {
        float dx = x - 3.7500000000000000e-01;
        return ((vec3( 2.3255546213942225e+02,  1.8349599099637384e+01,  1.7433813849989207e+01 ) * dx
               + vec3(-1.6448268217224928e+00, -2.9036166900100602e+00, -3.9744682262239461e+00 )) * dx
               + vec3(-5.1494864403514876e-02,  7.8968912818872505e-01, -8.7076475998425507e-01 )) * dx
               + vec3( 2.3090476190476199e-02,  6.4178571428571396e-01,  7.9126666666666701e-01 );
    } else if (x < 4.0625000000000000e-01) {
        float dx = x - 3.9062500000000000e-01;
        return ((vec3( 1.5126193200717549e+02,  2.0267550346934740e+01,  2.0857035135376179e+01 ) * dx
               + vec3( 9.2562104660629245e+00, -2.0434792322145579e+00, -3.1572582020057021e+00 )) * dx
               + vec3( 6.7433005039304356e-02,  7.1239075440396538e-01, -9.8219798542534331e-01 )) * dx
               + vec3( 2.2771428571428599e-02,  6.5348571428571400e-01,  7.7675714285714303e-01 );
    } else if (x < 4.2187500000000000e-01) {
        float dx = x - 4.0625000000000000e-01;
        return ((vec3( 1.0861181935568159e+02, -5.7969433444380156e+00,  3.9956456082908054e+00 ) * dx
               + vec3( 1.6346613528899276e+01, -1.0934378097019919e+00, -2.1795846800349437e+00 )) * dx
               + vec3( 4.6747712996058871e-01,  6.6337642562401933e-01, -1.0655861554572283e+00 )) * dx
               + vec3( 2.6661904761904800e-02,  6.6419523809523795e-01,  7.6071904761904796e-01 );
    } else if (x < 4.3750000000000000e-01) {
        float dx = x - 4.2187500000000000e-01;
        return ((vec3(-3.0484063800132168e+02,  1.4154965887634640e+01, -3.1353889969814710e+00 ) * dx
               + vec3( 2.1437792561196851e+01, -1.3651695289725239e+00, -1.9922887921463122e+00 )) * dx
               + vec3( 1.0578584751183406e+00,  6.2496068595722998e-01, -1.1307716784600605e+00 )) * dx
               + vec3( 3.8371428571428598e-02,  6.7427142857142897e-01,  7.4355238095238096e-01 );
    } else if (x < 4.5312500000000000e-01) {
        float dx = x - 4.3750000000000000e-01;
        return ((vec3( 1.9732370744832981e+01, -3.3873392535419122e+00, -5.1854420010455629e+00 ) * dx
               + vec3( 7.1483876548848961e+00, -7.0165550298965007e-01, -2.1392601513798186e+00 )) * dx
               + vec3( 1.5045175409946179e+00,  5.9266654483282100e-01, -1.1953271307026563e+00 )) * dx
               + vec3( 5.8971428571428598e-02,  6.8375714285714295e-01,  7.2538571428571397e-01 );
    } else if (x < 4.6875000000000000e-01) {
        float dx = x - 4.5312500000000000e-01;
        return ((vec3(-5.2460806882781675e+01, -6.0560887320505685e-01,  1.3890718905419471e+01 ) * dx
               + vec3( 8.0733425335489422e+00, -8.6043703049942721e-01, -2.3823277451788294e+00 )) * dx
               + vec3( 1.7423570751888966e+00,  5.6825884899705426e-01, -1.2659769415863851e+00 )) * dx
               + vec3( 8.4300000000000000e-02,  6.9283333333333297e-01,  7.0616666666666705e-01 );
    } else if (x < 4.8437500000000000e-01) {
        float dx = x - 4.6875000000000000e-01;
        return ((vec3( 1.0354971072183483e+01,  5.8097747460711062e+00, -5.4384621916749820e+00 ) * dx
               + vec3( 5.6142422109185510e+00, -8.8882494643091425e-01, -1.7312002964872917e+00 )) * dx
               + vec3( 1.9562255868212013e+00,  5.4092663060751767e-01, -1.3302508172374183e+00 )) * dx
               + vec3( 1.1329523809523800e-01,  7.0150000000000001e-01,  6.8585714285714305e-01 );
    } else if (x < 5.0000000000000000e-01) {
        float dx = x - 4.8437500000000000e-01;
        return ((vec3(-1.3925172644537971e+01, -8.9021377300786071e+00, -4.6199177582688593e+00 ) * dx
               + vec3( 6.0996314799271518e+00, -6.1649175520883115e-01, -1.9861282117220564e+00 )) * dx
               + vec3( 2.1392548632406654e+00,  5.1740605714439658e-01, -1.3883340751781894e+00 )) * dx
               + vec3( 1.4527142857142900e-01,  7.0975714285714298e-01,  6.6462857142857201e-01 );
    } else if (x < 5.1562500000000000e-01) {
        float dx = x - 5.0000000000000000e-01;
        return ((vec3( 3.1614367125520630e+01, -1.1395280968671647e+01,  2.1421523701702025e+01 ) * dx
               + vec3( 5.4468890122144344e+00, -1.0337794613062659e+00, -2.2026868566409092e+00 )) * dx
               + vec3( 2.3196692459303776e+00,  4.9162056938634824e-01, -1.4537843106213608e+00 )) * dx
               + vec3( 1.8013333333333301e-01,  7.1765714285714299e-01,  6.4243333333333297e-01 );
    } else if (x < 5.3125000000000000e-01) {
        float dx = x - 5.1562500000000000e-01;
        return ((vec3(-3.7634010143333590e+01,  2.0544616050328934e+00,  1.3219372364175872e+00 ) * dx
               + vec3( 6.9288124712232140e+00, -1.5679332567127493e+00, -1.1985529331236269e+00 )) * dx
               + vec3( 2.5130395816090907e+00,  4.5096880816730112e-01, -1.5069286823364316e+00 )) * dx
               + vec3( 2.1782857142857101e-01,  7.2504285714285699e-01,  6.1926190476190501e-01 );
    } else if (x < 5.4687500000000000e-01) {
        float dx = x - 5.3125000000000000e-01;
        return ((vec3( 1.2815768685879013e+01, -1.4298832118473902e+01,  3.9450879734146490e+01 ) * dx
               + vec3( 5.1647182457544520e+00, -1.4716303689768324e+00, -1.1365871251665525e+00 )) * dx
               + vec3( 2.7020009990618670e+00,  4.0347562651590141e-01, -1.5434152457472157e+00 )) * dx
               + vec3( 2.5864285714285701e-01,  7.3171428571428598e-01,  5.9542857142857097e-01 );
    } else if (x < 5.6250000000000000e-01) {
        float dx = x - 5.4687500000000000e-01;
        return ((vec3(-7.8540912219456771e+01, -1.8509114083431125e+01,  3.3113477160250433e+01 ) * dx
               + vec3( 5.7654574029050307e+00, -2.1418881245302965e+00,  7.1267286237156402e-01 )) * dx
               + vec3( 2.8727849935721714e+00,  3.4701440005485251e-01, -1.5500389061033872e+00 )) * dx
               + vec3( 3.0217142857142898e-01,  7.3760476190476199e-01,  5.7118571428571396e-01 );
    } else if (x < 5.7812500000000000e-01) {
        float dx = x - 5.6250000000000000e-01;
        return ((vec3(-5.8163891236508938e+01,  9.6920884524980497e+00,  3.0320583052976861e+01 ) * dx
               + vec3( 2.0838521426179946e+00, -3.0095028471911305e+00,  2.2648671042583031e+00 )) * dx
               + vec3( 2.9954304552209687e+00,  2.6652391612170523e-01, -1.5035148441247956e+00 )) * dx
               + vec3( 3.4816666666666701e-01,  7.4243333333333295e-01,  5.4726666666666701e-01 );
    } else if (x < 5.9375000000000000e-01) {
        float dx = x - 5.7812500000000000e-01;
        return ((vec3(-6.4543256167712116e+01, -2.8636353652780144e-01,  2.8905906284068501e+00 ) * dx
               + vec3(-6.4258025909336181e-01, -2.5551862009802844e+00,  3.6861444348665935e+00 )) * dx
               + vec3( 3.0179503284010409e+00,  1.7957564974402687e-01, -1.4105302888259692e+00 )) * dx
               + vec3( 3.9525714285714297e-01,  7.4590000000000001e-01,  5.2444285714285699e-01 );
    } else if (x < 6.0937500000000000e-01) {
        float dx = x - 5.9375000000000000e-01;
        return ((vec3(-2.4450284092939786e+01,  1.3922851408411924e+01, -1.6916850328844372e+01 ) * dx
               + vec3(-3.6680453919548675e+00, -2.5686094917550251e+00,  3.8216408705731646e+00 )) * dx
               + vec3( 2.9505968026034126e+00,  9.9516342045037676e-02, -1.2932211434284731e+00 )) * dx
               + vec3( 4.4200952380952402e-01,  7.4808095238095196e-01,  5.0331428571428605e-01 );
    } else if (x < 6.2500000000000000e-01) {
        float dx = x - 6.0937500000000000e-01;
        return ((vec3( 1.2547821111311350e+01,  1.5748329330961459e+01, -1.7611303598786566e+01 ) * dx
               + vec3(-4.8141524588114200e+00, -1.9159758319857161e+00,  3.0286635114085847e+00 )) * dx
               + vec3( 2.8180624611851890e+00,  2.9444696361588602e-02, -1.1861851374600081e+00 )) * dx
               + vec3( 4.8712380952380901e-01,  7.4906190476190504e-01,  4.8397619047619100e-01 );
    } else if (x < 6.4062500000000000e-01) {
        float dx = x - 6.2500000000000000e-01;
        return ((vec3( 9.2115329809656430e+00, -3.2661877796437579e+00, -1.2675733711774058e+00 ) * dx
               + vec3(-4.2259733442187004e+00, -1.1777728945968977e+00,  2.2031336552154643e+00 )) * dx
               + vec3( 2.6768104955128438e+00, -1.8895127491264742e-02, -1.1044383067315073e+00 )) * dx
               + vec3( 5.3002857142857096e-01,  7.4911428571428595e-01,  4.6611428571428598e-01 );
    } else if (x < 6.5625000000000000e-01) {
        float dx = x - 6.4062500000000000e-01;
        return ((vec3( 1.4269589821681299e+01,  7.3028598827757278e+00, -8.5260219639800940e+00 ) * dx
               + vec3(-3.7941827357359359e+00, -1.3308754467676989e+00,  2.1437161534415234e+00 )) * dx
               + vec3( 2.5514955567635522e+00, -5.8092757825086563e-02, -1.0365187784712420e+00 )) * dx
               + vec3( 5.7085714285714295e-01,  7.4851904761904797e-01,  4.4939047619047601e-01 );
    } else if (x < 6.7187500000000000e-01) {
        float dx = x - 6.5625000000000000e-01;
        return ((vec3( 8.6083934467238432e+00,  2.6914824850885094e-01, -1.7057138772896455e+01 ) * dx
               + vec3(-3.1252957128446250e+00, -9.8855388976258662e-01,  1.7440588738799565e+00 )) * dx
               + vec3( 2.4433787060044811e+00, -9.4333841208372265e-02, -9.7577229366934382e-01 )) * dx
               + vec3( 6.0985238095238103e-01,  7.4731428571428604e-01,  4.3368571428571401e-01 );
    } else if (x < 6.8750000000000000e-01) {
        float dx = x - 6.7187500000000000e-01;
        return ((vec3( 8.7188554392023345e+00,  1.7834947123447904e+01, -1.8886229447019101e+00 ) * dx
               + vec3(-2.7217772700294449e+00, -9.7593756561373424e-01,  9.4450549390043514e-01 )) * dx
               + vec3( 2.3520181906470738e+00, -1.2502902019862727e-01, -9.3376347542277516e-01 )) * dx
               + vec3( 6.4729999999999999e-01,  7.4560000000000004e-01,  4.1880000000000001e-01 );
    } else if (x < 7.0312500000000000e-01) {
        float dx = x - 6.8750000000000000e-01;
        return ((vec3( 8.9449847961700044e+00, -2.1676746266635202e+01, -4.0993789718798466e+00 ) * dx
               + vec3(-2.3130809213168355e+00, -1.3992441920211368e-01,  8.5597629336753311e-01 )) * dx
               + vec3( 2.2733485314072883e+00, -1.4246436371137491e-01, -9.0563094749671313e-01 )) * dx
               + vec3( 6.8341904761904804e-01,  7.4347619047619096e-01,  4.0443333333333298e-01 );
    } else if (x < 7.1875000000000000e-01) {
        float dx = x - 7.0312500000000000e-01;
        return ((vec3( 1.1674919661892304e+01,  2.3933066515154213e+01, -1.1673175453308831e+01 ) * dx
               + vec3(-1.8937847589963666e+00, -1.1560219004506387e+00,  6.6381790406066532e-01 )) * dx
               + vec3( 2.2076162551523946e+00, -1.6271352495594915e-01, -8.8188416316189755e-01 )) * dx
               + vec3( 7.1840952380952405e-01,  7.4113333333333298e-01,  3.9047619047618998e-01 );
    } else if (x < 7.3437500000000000e-01) {
        float dx = x - 7.1875000000000000e-01;
        return ((vec3(-4.4641682053710623e+00,  2.0910706819426692e+00,  4.6048045942407727e+00 ) * dx
               + vec3(-1.3465228998451648e+00, -3.4159407552784897e-02,  1.1663780468681384e-01 )) * dx
               + vec3( 2.1569864479829954e+00, -1.8131010789350266e-01, -8.6968954271271826e-01 )) * dx
               + vec3( 7.5248571428571398e-01,  7.3839999999999995e-01,  3.7681428571428599e-01 );
    } else if (x < 7.5000000000000000e-01) {
        float dx = x - 7.3437500000000000e-01;
        return ((vec3( 1.2423276968973711e+01, -6.0829492432479162e+00, -2.1725700066572116e+01 ) * dx
               + vec3(-1.5557807844719334e+00,  6.3859530663277708e-02,  3.3248802004185007e-01 )) * dx
               + vec3( 2.1116379529155407e+00, -1.8084604346990121e-01, -8.6267195170133282e-01 )) * dx
               + vec3( 7.8584285714285695e-01,  7.3556666666666704e-01,  3.6327142857142902e-01 );
    } else if (x < 7.6562500000000000e-01) {
        float dx = x - 7.5000000000000000e-01;
        return ((vec3( 3.4549460436900552e+00,  2.2240726291601970e+01, -7.5799471847609725e+00 ) * dx
               + vec3(-9.7343967655129060e-01, -2.2127871511396835e-01, -6.8590417057871789e-01 )) * dx
               + vec3( 2.0721188832120530e+00, -1.8330571822694325e-01, -8.6819407905347146e-01 )) * dx
               + vec3( 8.1850476190476196e-01,  7.3273333333333301e-01,  3.4979047619047599e-01 );
    } else if (x < 7.8125000000000000e-01) {
        float dx = x - 7.6562500000000000e-01;
        return ((vec3( 8.7094721894791203e+00,  1.3239510743088688e+01, -2.2852796908624047e+01 ) * dx
               + vec3(-8.1148908075331927e-01,  8.2125532980487381e-01, -1.0412141948643885e+00 )) * dx
               + vec3( 2.0442293713791684e+00, -1.7393108362239784e-01, -8.9518030351351996e-01 )) * dx
               + vec3( 8.5065714285714300e-01,  7.2989999999999999e-01,  3.3602857142857101e-01 );
    } else if (x < 7.9687500000000000e-01) {
        float dx = x - 7.8125000000000000e-01;
        return ((vec3(-1.2078434801289291e+01,  4.3390183117236198e+01, -3.9570693752303733e+01 ) * dx
               + vec3(-4.0323257187148548e-01,  1.4418573958871561e+00, -2.1124390499561407e+00 )) * dx
               + vec3( 2.0252493455569058e+00, -1.3856994728345987e-01, -9.4445613546384066e-01 )) * dx
               + vec3( 8.8243333333333296e-01,  7.2743333333333304e-01,  3.2169999999999999e-01 );
    } else if (x < 8.1250000000000000e-01) {
        float dx = x - 7.9687500000000000e-01;
        return ((vec3(-1.2824532984374384e+01,  1.1653781393088177e+02, -1.1096774236821523e+02 ) * dx
               + vec3(-9.6940920318192092e-01,  3.4757722295076028e+00, -3.9673153195953783e+00 )) * dx
               + vec3( 2.0038018178216963e+00, -6.1731984386666772e-02, -1.0394522974880831e+00 )) * dx
               + vec3( 9.1393333333333304e-01,  7.2578571428571403e-01,  3.0627619047619098e-01 );
    } else if (x < 8.2812500000000000e-01) {
        float dx = x - 8.1250000000000000e-01;
        return ((vec3(-3.5855044278532131e+02,  2.7064903734930277e+02, -8.0792089155266083e+01 ) * dx
               + vec3(-1.5705591868244702e+00,  8.9384822575176859e+00, -9.1689282431054675e+00 )) * dx
               + vec3( 1.9641148117278464e+00,  1.3224074197310332e-01, -1.2447061031552840e+00 )) * dx
               + vec3( 9.4495714285714305e-01,  7.2611428571428605e-01,  2.8864285714285698e-01 );
    } else if (x < 8.4375000000000000e-01) {
        float dx = x - 8.2812500000000000e-01;
        return ((vec3(-3.8174017206443654e+02, -1.9549693475620506e+02,  4.4911575613188438e+02 ) * dx
               + vec3(-1.8377611192386407e+01,  2.1625155883266252e+01, -1.2956057422258565e+01 )) * dx
               + vec3( 1.6524246495526764e+00,  6.0979758792285232e-01, -1.5904090041765968e+00 )) * dx
               + vec3( 9.7389523809523804e-01,  7.3139523809523799e-01,  2.6664761904761902e-01 );
    } else if (x < 8.5937500000000000e-01) {
        float dx = x - 8.4375000000000000e-01;
        return ((vec3( 4.3248438818547703e+02, -2.7134838403902307e+02,  3.3204036056432756e+01 ) * dx
               + vec3(-3.6271681757906869e+01,  1.2461237066569140e+01,  8.0962436464235150e+00 )) * dx
               + vec3( 7.9852944720434427e-01,  1.1423974777640304e+00, -1.6663435944240195e+00 )) * dx
               + vec3( 9.9377142857142897e-01,  7.4545714285714304e-01,  2.4034761904761900e-01 );
    } else if (x < 8.7500000000000000e-01) {
        float dx = x - 8.5937500000000000e-01;
        return ((vec3( 1.7847934313241271e+02, -6.1117386114828536e+00, -1.0882439559595376e+02 ) * dx
               + vec3(-1.5998976061712632e+01, -2.5821843526006538e-01,  9.6526828365688004e+00 )) * dx
               + vec3(-1.8199581227210410e-02,  1.3330696438782346e+00, -1.3890166181272647e+00 )) * dx
               + vec3( 9.9904285714285701e-01,  7.6531428571428595e-01,  2.1641428571428600e-01 );
    } else if (x < 8.9062500000000000e-01) {
        float dx = x - 8.7500000000000000e-01;
        return ((vec3( 1.0065469642774150e+02,  1.1181852770679304e+01, -4.2302948910418884e+01 ) * dx
               + vec3(-7.6327568523807861e+00, -5.4470618267332416e-01,  4.5515392930084682e+00 )) * dx
               + vec3(-3.8744540800992006e-01,  1.3205239467230254e+00, -1.1670756473526198e+00 )) * dx
               + vec3( 9.9553333333333305e-01,  7.8605714285714301e-01,  1.9665238095238100e-01 );
    } else if (x < 9.0625000000000000e-01) {
        float dx = x - 8.9062500000000000e-01;
        return ((vec3( 5.1792385442186948e+01,  1.3813127528788970e+01, -4.7771351619749993e+01 ) * dx
               + vec3(-2.9145679573304033e+00, -2.0556834047731776e-02,  2.5685885628325829e+00 )) * dx
               + vec3(-5.5224735816165738e-01,  1.3116917120867588e+00, -1.0558236496051034e+00 )) * dx
               + vec3( 9.8799999999999999e-01,  8.0659999999999998e-01,  1.7936666666666701e-01 );
    } else if (x < 9.2187500000000000e-01) {
        float dx = x - 9.0625000000000000e-01;
        return ((vec3( 1.1035785704157649e+02,  5.2154589495154021e+01, -3.9990387467675163e+01 ) * dx
               + vec3(-4.8679988972789023e-01,  6.2693351886425119e-01,  3.2930645565680206e-01 )) * dx
               + vec3(-6.0539373077194325e-01,  1.3211663477870170e+00, -1.0105440399412067e+00 )) * dx
               + vec3( 9.7885714285714298e-01,  8.2714285714285696e-01,  1.6331428571428599e-01 );
    } else if (x < 9.3750000000000000e-01) {
        float dx = x - 9.2187500000000000e-01;
        return ((vec3( 4.6043843534396274e+01,  2.0987943062129727e+01, -2.3203479461840441e+01 ) * dx
               + vec3( 4.6862246590960082e+00,  3.0716799014495959e+00, -1.5452429568904713e+00 )) * dx
               + vec3(-5.3977771875056635e-01,  1.3789571824794209e+00, -1.0295430477729828e+00 )) * dx
               + vec3( 9.6970000000000001e-01,  8.4813809523809502e-01,  1.4745238095238100e-01 );
    } else if (x < 9.5312500000000000e-01) {
        float dx = x - 9.3750000000000000e-01;
        return ((vec3( 6.1233625963980650e+01,  2.8669866827404956e+01,  2.4201791029260814e+01 ) * dx
               + vec3( 6.8445298247708335e+00,  4.0554897324869268e+00, -2.6329060566642419e+00 )) * dx
               + vec3(-3.5960967994014698e-01,  1.4903192080096790e+00, -1.0948266261097752e+00 )) * dx
               + vec3( 9.6258571428571404e-01,  8.7051428571428602e-01,  1.3089999999999999e-01 );
    } else if (x < 9.6875000000000000e-01) {
        float dx = x - 9.5312500000000000e-01;
        return ((vec3( 4.1070719275903762e+01,  5.3910277236601019e+00,  2.0019172487757277e+01 ) * dx
               + vec3( 9.7148560418324266e+00,  5.3993897400215340e+00, -1.4984471021676413e+00 )) * dx
               + vec3(-1.0086927577447102e-01,  1.6380516997676238e+00, -1.1593790192165234e+00 )) * dx
               + vec3( 9.5887142857142904e-01,  8.9490000000000003e-01,  1.1324285714285701e-01 );
    } else if (x < 9.8437500000000000e-01) {
        float dx = x - 9.6875000000000000e-01;
        return ((vec3(-5.3250445924665847e+01, -1.6529749150400146e+01, -1.4422423336140781e+02 ) * dx
               + vec3( 1.1640046007890415e+01,  5.6520941645681013e+00, -5.6004839180401900e-01 )) * dx
               + vec3( 2.3280106875244833e-01,  1.8107311357768368e+00, -1.1915430113098306e+00 )) * dx
               + vec3( 9.5982380952380997e-01,  9.2183333333333295e-01,  9.4838095238095305e-02 );
    } else if (x < 1.0000000000000000e+00) {
        float dx = x - 9.8437500000000000e-01;
        return ((vec3(-1.9507053557699635e+02, -1.0404825969371934e+02,  1.5617193238656020e+02 ) * dx
               + vec3( 9.1439313551717039e+00,  4.8772621731430945e+00, -7.3205593306200099e+00 )) * dx
               + vec3( 5.5755071505029385e-01,  1.9752523285535741e+00, -1.3146775069727061e+00 )) * dx
               + vec3( 9.6609999999999996e-01,  9.5144285714285703e-01,  7.5533333333333300e-02 );
    } else {
        float dx = x - 1.0000000000000000e+00;
        return ((vec3( 0.0000000000000000e+00,  3.4202936336155174e+00,  3.0625241907655076e+00 ) * dx
               + vec3( 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00 )) * dx
               + vec3( 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00 )) * dx
               + vec3( 9.7629999999999995e-01,  9.8309999999999997e-01,  5.3800000000000001e-02 );
    }
}*/

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
    float len = std::max( abs(min), abs(max));
    if(Ks<0)
    {
        return (len-Ks)/(2*len);
    }
    else
    {
        return (len+Ks)/(2*len);
    }
}

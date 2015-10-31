## Plotting

# GLPlot is already imported
import ModernGL
import GLAbstraction
import GLFW

#export surf

function colorf(ext;colormap::Symbol=:hottocold)
    if colormap == :hottocold
        hottocold(ext[2],ext[1])
    elseif colormap == :cubehelix
        cubehelix(ext[2],ext[1])
    end
end

#
# This colormap is based on the CubeHelix colormap described in
# D. A. Green, "A colour scheme for the display of astronomical intensity images," Bull. Astr. Soc. India (2011) 39, 289–295.
# "This – unlike many currently available schemes – is designed to be monotonically increasing in terms of its perceived brightness."
#
cubehelix(ma,mi) = "vec4(   (xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*"))*( 1.0+(1.0-(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))/2.0*(-0.14861*cos(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )+1.78277*sin(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )) )   ,    (xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*"))*( 1.0+(1.0-(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))/2.0*(-0.29227*cos(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )-0.90649*sin(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )) )       ,   (xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*"))*( 1.0+(1.0-(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))/2.0*(1.97294*cos(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )) )    ,0.5 );"

#
# This hot-to-cold colormap is described in
# http://paulbourke.net/texture_colour/colourspace/
#
hottocold(ma,mi) = "xyz.z > (("*string(mi)*")+0.75*(("*string(ma)*")-("*string(mi)*"))) ? vec4(1.0, 1.0 + 4.0*(("*string(mi)*")+0.75*(("*string(ma)*")-("*string(mi)*")) - xyz.z)/(("*string(ma)*")-("*string(mi)*"))  ,0.0,0.5) :  xyz.z > (("*string(mi)*")+0.5*(("*string(ma)*")-("*string(mi)*"))) ? vec4( 4.0*(xyz.z-("*string(mi)*")-0.5*(("*string(ma)*")-("*string(mi)*")))/(("*string(ma)*")-("*string(mi)*")) ,1.0,0.0,0.5) : xyz.z > (("*string(mi)*")+0.25*(("*string(ma)*")-("*string(mi)*"))) ? vec4(0.0,1.0,1.0 + 4.0*(("*string(mi)*")+0.25*(("*string(ma)*")-("*string(mi)*"))-xyz.z)/(("*string(ma)*")-("*string(mi)*")),0.5) : vec4(0.0,4.0*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")),1.0,0.5);"

function glupdatewindow(obj,window)
    ModernGL.glClear(ModernGL.GL_COLOR_BUFFER_BIT | ModernGL.GL_DEPTH_BUFFER_BIT)
    GLAbstraction.render(obj)
    # GLFW.SwapBuffers(window.glfwWindow)
    GLFW.PollEvents()
    yield()
    obj,window
end


## Vector routines
function glsurfupdate(vals::Matrix,obj,window)##obj should be type RenderObject, window should be type Screen
    zvalues = obj.uniforms[:z]
#    colrs=obj.uniforms[:color]
    zvalues[:,:]=float32(vals)
#    color     = map(colorf,vals)
#    GLAbstraction.update!(zvalues, map(GLAbstraction.Vec1,vals)) # now you can simply update the gpu object, which should be very efficient
#     GLAbstraction.update!(colrs,color)

    glupdatewindow(obj,window)
end

function glsurf(vals::Matrix;colormap::Symbol=:hottocold)
    window = GLPlot.createdisplay(w=1000,h=1000,eyeposition=GLAbstraction.Vec3(1.,1.,1.), lookat=GLAbstraction.Vec3(0.,0.,0.),async=true)
    ModernGL.glClearColor(1,1,1,0)
    obj     = GLPlot.glplot(map(GLAbstraction.Vec1,vals) , primitive=GLPlot.SURFACE(), color=colorf(extrema(vals);colormap=colormap))


    glupdatewindow(obj,window)
end


function glsurf(xx::Array,yy::Array,vals::Matrix;colormap::Symbol=:hottocold)
    window = GLPlot.createdisplay(w=1000,h=1000,eyeposition=GLAbstraction.Vec3(1.,1.,1.), lookat=GLAbstraction.Vec3(0.,0.,0.),async=true)
    ModernGL.glClearColor(1,1,1,0)
    obj     = GLPlot.glplot(map(GLAbstraction.Vec1,vals) , xrange=float32(xx),yrange=float32(yy),primitive=GLPlot.SURFACE(), color=colorf(extrema(vals);colormap=colormap))


    glupdatewindow(obj,window)
end










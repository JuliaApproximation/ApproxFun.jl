## Plotting


#export surf

#
# This colormap is based on the CubeHelix colormap described in
# D. A. Green, "A colour scheme for the display of astronomical intensity images," Bull. Astr. Soc. India (2011) 39, 289–295.
# "This – unlike many currently available schemes – is designed to be monotonically increasing in terms of its perceived brightness."
#
colorf(ma,mi) = "vec4(   (xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*"))*( 1.0+(1.0-(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))/2.0*(-0.14861*cos(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )+1.78277*sin(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )) )   ,    (xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*"))*( 1.0+(1.0-(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))/2.0*(-0.29227*cos(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )-0.90649*sin(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )) )       ,   (xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*"))*( 1.0+(1.0-(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))/2.0*(1.97294*cos(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )+0.0*sin(  6.28*(0.5/3.0+1.0-1.5*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")))  )) )    ,0.45 +0.0*(xyz.z-("*string(mi)*"))/(("*string(ma)*")-("*string(mi)*")));"

function glupdatewindow(obj,window)
    require("GLPlot")
    ModernGL=Main.ModernGL
    GLAbstraction=Main.GLAbstraction

     ModernGL.glClear(ModernGL.GL_COLOR_BUFFER_BIT | ModernGL.GL_DEPTH_BUFFER_BIT)
    GLAbstraction.render(obj)
#    Main.GLFW.SwapBuffers(window.glfwWindow)
    Main.GLFW.PollEvents()
    yield()
    obj,window
end


## Vector routines
function glsurfupdate(vals::Matrix,obj,window)##obj should be type RenderObject, window should be type Screen
    require("GLPlot")
    GLAbstraction=Main.GLAbstraction

    zvalues = obj.uniforms[:z]
#    colrs=obj.uniforms[:color]
    zvalues[:,:]=float32(vals)
#    color     = map(colorf,vals)
#    GLAbstraction.update!(zvalues, map(GLAbstraction.Vec1,vals)) # now you can simply update the gpu object, which should be very efficient
#     GLAbstraction.update!(colrs,color)

    glupdatewindow(obj,window)
end

function glsurf(vals::Matrix)
    require("GLPlot")
    GLAbstraction=Main.GLAbstraction
    ModernGL=Main.ModernGL
    GLPlot=Main.GLPlot

    window = GLPlot.createdisplay(w=1000,h=1000,eyeposition=GLAbstraction.Vec3(1.,1.,1.), lookat=GLAbstraction.Vec3(0.,0.,0.),async=true)
    mi,ma=extrema(vals)
    ModernGL.glClearColor(1,1,1,0)
    obj     = GLPlot.glplot(map(GLAbstraction.Vec1,vals) , primitive=GLPlot.SURFACE(), color=colorf(ma,mi))


    glupdatewindow(obj,window)
end


function glsurf(xx::Array,yy::Array,vals::Matrix)
    require("GLPlot")
    GLAbstraction=Main.GLAbstraction
    ModernGL=Main.ModernGL
    GLPlot=Main.GLPlot

    window = GLPlot.createdisplay(w=1000,h=1000,eyeposition=GLAbstraction.Vec3(1.,1.,1.), lookat=GLAbstraction.Vec3(0.,0.,0.),async=true)
    mi,ma=extrema(vals)
    ModernGL.glClearColor(1,1,1,0)
    obj     = GLPlot.glplot(map(GLAbstraction.Vec1,vals) , xrange=float32(xx),yrange=float32(yy),primitive=GLPlot.SURFACE(), color=colorf(ma,mi))


    glupdatewindow(obj,window)
end










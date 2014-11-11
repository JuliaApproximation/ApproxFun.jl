## Plotting


#export surf


colorf(m)="xyz.z>0 ? vec4(.2,.4,0.5+"*string(m)*"*xyz.z,0.98) : vec4(.2,.4-"*string(1.2m)*"*xyz.z,0.5,0.98);"


function glupdatewindow(obj,window)
    ModernGL=Main.ModernGL
    GLAbstraction=Main.GLAbstraction
        
     ModernGL.glClear(ModernGL.GL_COLOR_BUFFER_BIT | ModernGL.GL_DEPTH_BUFFER_BIT)
    GLAbstraction.render(obj)
    Main.GLFW.SwapBuffers(window.glfwWindow)     
    Main.GLFW.PollEvents()
    yield()    
    obj,window   
end


## Vector routines
function glsurfupdate(vals::Matrix,obj,window)##obj should be type RenderObject, window should be type Screen
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
    GLAbstraction=Main.GLAbstraction
    ModernGL=Main.ModernGL
    GLPlot=Main.GLPlot
    
    window = GLPlot.createdisplay(w=1000,h=1000,eyeposition=GLAbstraction.Vec3(1.,1.,1.), lookat=GLAbstraction.Vec3(0.,0.,0.),async=true)
    m=0.5./maximum(abs(vals))
    ModernGL.glClearColor(1,1,1,0)
    obj     = GLPlot.glplot(map(GLAbstraction.Vec1,vals) , primitive=GLPlot.SURFACE(), color=colorf(m))


    glupdatewindow(obj,window)
end


function glsurf(xx::Array,yy::Array,vals::Matrix)
    GLAbstraction=Main.GLAbstraction
    ModernGL=Main.ModernGL
    GLPlot=Main.GLPlot
    
    window = GLPlot.createdisplay(w=1000,h=1000,eyeposition=GLAbstraction.Vec3(1.,1.,1.), lookat=GLAbstraction.Vec3(0.,0.,0.),async=true)
    m=0.5./maximum(abs(vals))  
    ModernGL.glClearColor(1,1,1,0)
    obj     = GLPlot.glplot(map(GLAbstraction.Vec1,vals) , xrange=float32(xx),yrange=float32(yy),primitive=GLPlot.SURFACE(), color=colorf(m))


    glupdatewindow(obj,window)
end










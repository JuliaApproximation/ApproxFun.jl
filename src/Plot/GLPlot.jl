## Plotting


export surf


colorf(x)=x>0?Main.GLAbstraction.Vec4(.1,.1,0.5+3x,.1):Main.GLAbstraction.Vec4(.1,.1-3x,0.5,.1)


function glupdatewindow(obj,window)
    ModernGL=Main.ModernGL
    GLAbstraction=Main.GLAbstraction
        
     ModernGL.glClear(ModernGL.GL_COLOR_BUFFER_BIT | ModernGL.GL_DEPTH_BUFFER_BIT)
    GLAbstraction.render(obj)
    Main.GLFW.SwapBuffers(window.glfwWindow)     
    Main.GLFW.PollEvents()
    obj,window   
end


## Vector routines
function surf(vals::Matrix,obj,window)##obj should be type RenderObject, window should be type Screen
    GLAbstraction=Main.GLAbstraction
    
    zvalues = obj.uniforms[:z] 
    colrs=obj.uniforms[:color]    
    
    color     = map(colorf,vals)    
    GLAbstraction.update!(zvalues, map(GLAbstraction.Vec1,vals)) # now you can simply update the gpu object, which should be very efficient
     GLAbstraction.update!(colrs,color)

    glupdatewindow(obj,window)
end

function surf(vals::Matrix)
    GLAbstraction=Main.GLAbstraction
    ModernGL=Main.ModernGL
    GLPlot=Main.GLPlot
    
    window = GLPlot.createdisplay(w=1000,h=1000,eyeposition=GLAbstraction.Vec3(1.,1.,1.), lookat=GLAbstraction.Vec3(0.,0.,0.))
    
    ModernGL.glClearColor(1,1,1,0)
    obj     = GLPlot.glplot(map(GLAbstraction.Vec1,vals) , primitive=GLPlot.SURFACE(), color=map(colorf,vals))


    glupdatewindow(obj,window)
end




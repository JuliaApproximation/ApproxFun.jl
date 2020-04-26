using ApproxFun, Test, DualNumbers
import ApproxFun: eigs

@testset "Extras" begin
    @testset "Dual numbers" begin
        @test dual(1.5,1) ∈  Segment(dual(1.0,1),dual(2.0))

        f=Fun(exp,Segment(dual(1.0,1),dual(2.0)),20)
        @test Fun(h->Fun(exp,Segment(1.0+h,2.0)).coefficients[1],0..1)'(0.) ≈ DualNumbers.epsilon(f.coefficients[1])


        ud=let d= Segment(dual(0.0,1.0),1.0)
            B = ldirichlet(d)
            D = Derivative(d)
            a = Fun(exp,d)
            u = [B;D+dual(1.0,4.0)*a] \ [dual(1.0,2.0),0.0]
            u(0.5)
        end

        u0=let d= Segment(0.0,1.0)
            B = ldirichlet(d)
            D = Derivative(d)
            a = Fun(exp,d)
            u = [B;D+a] \ [1.0,0.0]
            u(0.5)
        end
        h=0.00001
        uh=let d= Segment(h,1.0)
            B = ldirichlet(d)
            D = Derivative(d)
            a = Fun(exp,d)
            u = [B;D+(1+4h)*a] \ [1.0+2h,0.0]
            u(0.5)
        end

        @test absdual(ud - dual(u0,(uh-u0)/h)) ≤ h

        let d=Segment(0.0,1.0)
            B = ldirichlet(d)
            D = Derivative(d)
            a = Fun(exp,d)
            u = [B;D+a] \ [dual(1.0,2.0),0.0]
            ur = [B;D+a] \ [1.0,0.0]
            ud = [B;D+a] \ [2.0,0.0]
            @test absdual(u(0.5)  - dual(ur(0.5),ud(0.5))) < 10eps()
        end
    end

    @testset "Eig test #336" begin
        d = 0..π
        A=Derivative(d)^2
        λ=eigvals(Dirichlet(d),A,100)
        @test sort(λ)[end-5:end] ≈ -(-6:-1).^2
        λ=eigvals(Dirichlet(), A, 100)
        @test sort(λ)[end-5:end] ≈ -(-6:-1).^2

        F = x->x^8
        d = Interval(0.0,1.0)
        f = Fun(F,d)
        ginf = Fun(x->exp(-x),d)
        gp = ginf'
        Af = Fun(x->x+f(x),d)
        transport_ = Fun(x-> x - 1,d)
        damping = Fun(x-> 1 - f(x),d)
        A = transport_*Derivative(d) + damping
        P = -DefiniteIntegral(Chebyshev(d))[LowRankFun((x,y)->gp(x)*(y+f(y)),d^2)];
        λ,V = eigs(A,100)
        @test norm(sort(real(filter(x->isreal(x),λ)))[1:5]-(0:4)) ≤ 100000eps()

        λ,V = eigs(A+P,100)
        @test sort(real(filter(x->isreal(x),λ)))[5] ≈ 3.93759261234502 atol=1E-3
    end

    @testset "findmin/max" begin
        f = Fun(x -> exp(0.25x) + sin(x) + 0.5cos(10x), -4..4)
        @test [findmax(f)...] ≈ [3.0531164509549584, 1.886754631165656]
        @test [findmin(f)...] ≈ [-0.825047261209411, -1.5741041425422948]
    end

    @testset "Newton iteration" begin
        x=Fun(Chebyshev(0..1))

        # Test two coupled systems of nonlinear ODEs with nonlinear BCs
        u1 = 0.5*one(x)
        u2 = 0.5*one(x)

        N_dep = (u1,u2) -> [
                u1(0)*u1(0)*u1(0) - 0.125;
                u2(0) - 1;
                u1' - u1*u1;  
                u2' + u1*u1*u2;
            ] 
        u1,u2 = newton(N_dep, [u1,u2])

        u1_exact = -1 / (x - 2)
        u2_exact = exp(1 / (x - 2) + 1/2)
        @test norm(u2 - u2_exact) < 1e-15
        @test norm(u1 - u1_exact) < 1e-15

        # Test two independent systems of nonlinear ODEs with nonlinear BCs
        u1 = 0.5*one(x)
        u2 = 0.5*one(x)

        N_ind = (u1,u2) -> [
                u1(0)*u1(0)*u1(0) - 0.125;
                u2(0)*u2(0)*u2(0) + 0.125;
                u1' - u1*u1;  
                u2' - u2*u2;
            ] 

        # note takes a few more iterations to converge to accuracy
        u1,u2 = newton(N_ind, [u1,u2], maxiterations=25) 

        u1_exact = -1 / (x - 2)
        u2_exact = -1 / (x + 2)
        @test norm(u2 - u2_exact) < 1e-15
        @test norm(u1 - u1_exact) < 1e-15

        # Compare to original newton for one equation and verify all solutions correct
        N_single = u -> [u(0)*u(0)*u(0) - 0.125;
                        u' - u*u]

        u = 0.5*one(x)
        u = newton(N_single, [u])
        @test norm(u - u1) < 1e-15
    end

    @testset "cumsum (#701)" begin
        d = ApproxFun.DualFun(Fun())
        @test cumsum(d).J*Fun() ≈ Fun(x -> (x^2-1)/2)
        @test cumsum(d*d).J*Fun(1) ≈ cumsum(2d).f ≈ Fun(x -> x^2-1)
    end
end

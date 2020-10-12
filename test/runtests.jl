
using SafeTestsets



@safetestset "Probability Objects" begin
    using InformationGeometry, Test, LinearAlgebra, Distributions

    DS = DataSet([0,0.5,1],[1.,3.,7.],[1.2,2.,0.6])
    model(x,p) = p[1] .* x .+ p[2]
    DM = DataModel(DS,model)
    XYPlane = Plane([0,0,0],[1,0,0],[0,1,0])
    x = [0.1,0.5]; R = rand(2)
    p = rand(2)

    @test IsLinear(DM)
    Dist = product_distribution([Normal(ydata(DM)[i],sigma(DM)[i]) for i in 1:length(ydata(DM))])
    @test loglikelihood(DM,p) ≈ logpdf(Dist,EmbeddingMap(DM,p))
    @test Score(DM,p) ≈ transpose(EmbeddingMatrix(DM,p)) * gradlogpdf(Dist,EmbeddingMap(DM,p))
    @test FisherMetric(DM,p) ≈ transpose(EmbeddingMatrix(DM,p)) * inv(cov(Dist)) * EmbeddingMatrix(DM,p)

    # Test AD vs manual derivative
    @test sum(abs.(InformationGeometry.AutoScore(DM,p) .- Score(DM,p))) < 2e-13
    @test sum(abs.(InformationGeometry.AutoMetric(DM,p) .- FisherMetric(DM,p))) < 2e-9

    # Do these tests in higher dimensions, check that OrthVF(PL) IsOnPlane....
    # @test OrthVF(DM,XYPlane,p) == OrthVF(DM,p)
    @test dot(OrthVF(DM,p),Score(DM,p)) < 2e-15
    @test sum(abs.(FindMLE(DM) .- [6.1213483146067,0.8382022471910])) < 1e-10
    # ALSO DO NONLINEAR MODEL!

    sols = MultipleConfidenceRegions(DM,1:2)
    @test StructurallyIdentifiable(DM,sols[1])[1] == true
    @test size(SaveConfidence(sols,100)) == (100,4)
    @test size(SaveGeodesics(sols,100)) == (100,2)
end


@safetestset "Inputting Dataset of various shapes" begin
    using InformationGeometry, Test, LinearAlgebra

    # Covariance Matrices, multidimensional Y, multidimensional X
    @test xdata(DataSet([[0,1.],[2,3],[5,8]],[1,5,6],[0.7,1.1,2.2])) == [[0,1.],[2,3],[5,8]]
    @test ydata(DataSet([1,5,6],[[0,1.],[2,3],[5,8]],[0.7,1.1,2.2])) == [[0,1.],[2,3],[5,8]]
    
end


# Test DataSets of different ydim and xdim and non-linear models.


@safetestset "Kullback-Leibler Divergences" begin
    using InformationGeometry, Test, LinearAlgebra, Distributions

    # Analytical divergences via types defined in Distributions.jl
    @test KullbackLeibler(MvNormal([1,2,3.],diagm([2,4,5.])),MvNormal([1,-2,-3.],diagm([1,5,2.]))) ≈ 11.056852819440055
    @test KullbackLeibler(Normal(1,2),Normal(5,1)) ≈ 8.806852819440055
    @test KullbackLeibler(Cauchy(1,3),Cauchy(5,2)) ≈ 0.5355182363563621
    @test KullbackLeibler(Exponential(5),Exponential(1)) ≈ 2.3905620875658995
    @test KullbackLeibler(Weibull(12,2),Weibull(3,4)) ≈ 2.146124463755512
    @test KullbackLeibler(Gamma(5,15),Gamma(20,1)) ≈ 29.409061308330323


    # Numerically calculated for via arbitrary types defined in Distributions.jl
    # ALSO ADD TESTS FOR DISCRETE DISTRIBUTIONS, DISTRIBUTIONS WITH LIMITED DOMAIN
    @test abs(KullbackLeibler(Cauchy(1,2),Normal(5,1),HyperCube([[-20,20]]),Carlo=false) - 16.77645704773449) < 1e-9
    @test abs(KullbackLeibler(Cauchy(1,2),Normal(5,1),HyperCube([[-20,20]]),Carlo=true,N=Int(3e6)) - 16.7764) < 5e-2
    @test abs(KullbackLeibler(MvTDist(1,[3,2,1.],diagm([1.,2.,3.])),MvNormal([1,2,3.],diagm([2,4,5.])),HyperCube(collect([-10,10.] for i in 1:3)),N=Int(3e6)) - 1.6559288) < 3e-1

    # Via any positive (hopefully normalized) functions
    @test abs(KullbackLeibler(x->pdf.(Normal(1,3),x),y->pdf.(Normal(5,2),y),Carlo=true,N=Int(3e6)) - 2.2195303) < 2e-2
    @test abs(KullbackLeibler(x->pdf.(Normal(1,3),x),y->pdf.(Normal(5,2),y),Carlo=false) - 2.21953032578115) < 1e-9
    @test abs(KullbackLeibler(x->pdf(MvNormal([1,2,3.],diagm([1,2,1.5])),x),y->pdf(MvNormal([1,-2,-3.],diagm([2,1.5,1.])),y),HyperCube(collect([-7,7.] for i in 1:3)),N=Int(3e6)) - 23.4) < 5e-1

end

@safetestset "Differential Geometry" begin
    using InformationGeometry, Test

    function S2metric(θ,ϕ)
        metric = zeros(typeof(ϕ),2,2);    metric[1,1] = 1.;    metric[2,2] = sin(θ)^2
        metric
    end
    S2metric(p::Vector) = S2metric(p...)

    function S2Christoffel(θ,ϕ)
        Symbol = zeros(typeof(ϕ),2,2,2);    Symbol[1,2,2] = -sin(θ)*cos(θ)
        Symbol[2,1,2] = cot(θ);    Symbol[2,2,1] = cot(θ)
        Symbol
    end
    S2Christoffel(p::Vector) = S2Christoffel(p...)
    ConstMetric(x) = [1. 0.; 0. 1.]

    # Test Numeric Christoffel Symbols, Riemann and Ricci tensors, Ricci Scalar
    # Test WITH AND WITHOUT BIGFLOAT
    x = rand(2)
    @test sum(abs.(ChristoffelSymbol(S2metric,x) .- S2Christoffel(x))) < 5e-10
    @test sum(abs.(ChristoffelSymbol(S2metric,BigFloat.(x)) .- S2Christoffel(BigFloat.(x)))) < 1e-40

    @test abs(RicciScalar(S2metric,rand(2)) - 2) < 5e-4
    @test abs(RicciScalar(S2metric,rand(BigFloat,2)) - 2) < 2e-22

    @test abs(GeodesicDistance(ConstMetric,[0,0],[1,1]) - sqrt(2)) < 1e-13
    @test abs(GeodesicDistance(S2metric,[pi/4,1],[3pi/4,1]) - pi/2) < 1e-11
    @test abs(GeodesicDistance(S2metric,[pi/2,0],[pi/2,pi/2]) - pi/2) < 3e-10

    # DS = DataSet([0,0.5,1],[1.,3.,7.],[1.2,2.,0.6])
    # model(x,p) = p[1].^3 .*x .+ p[2].^3
    # DM = DataModel(DS,model);   MLE = FindMLE(DM)
    # geo = GeodesicBetween(DM,MLE,MLE + rand(2),tol=1e-11)
    # @test sum(abs.(MLE .- [1.829289173660125,0.942865200406147])) < 1e-7
    # @test abs(InformationGeometry.ParamVol(geo) * InformationGeometry.GeodesicEnergy(DM,geo) - GeodesicLength(DM,geo)^2) < 1e-8
end


@safetestset "Numerical Helper Functions" begin
    using InformationGeometry, Test

    # Test integration, differentiation, Monte Carlo, GeodesicLength
    # TEST WITH AND WITHOUT BIGFLOAT
    @test abs(InformationGeometry.MonteCarloArea(x->((x[1]^2 + x[2]^2) < 1), HyperCube([[-1,1],[-1,1]])) - pi) < 1.5e-3
    @test abs(Integrate1D(cos,[0,pi/2]) .- 1) < 1e-13
    z = 3rand()
    @test abs(Integrate1D(x->2/sqrt(pi) * exp(-x^2),[0,z/sqrt(2)]) - ConfVol(z)) < 1e-12
    @test abs(LineSearch(x->(x < BigFloat(pi))) - pi) < 1e-14
    @test abs(CubeVol(TranslateCube(HyperCube([[0,1],[0,pi],[-sqrt(2),0]]),rand(3))) - sqrt(2)*pi) < 3e-15
end

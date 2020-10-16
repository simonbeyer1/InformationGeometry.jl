var documenterSearchIndex = {"docs":
[{"location":"confidence-regions/#Confidence-Regions","page":"Confidence Regions","title":"Confidence Regions","text":"","category":"section"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"Once a DataModel object has been defined, it can subsequently be used to compute various quantities as follows:","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"loglikelihood(::DataModel,::Vector{Float64})\nMLE(::DataModel)\nLogLikeMLE(::DataModel)","category":"page"},{"location":"confidence-regions/#StatsBase.loglikelihood-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"StatsBase.loglikelihood","text":"loglikelihood(DM::DataModel, θ::Vector) -> Real\n\nCalculates the logarithm of the likelihood L, i.e. ell(mathrmdata    theta) coloneqq mathrmln big( L(mathrmdata    theta) big) given a DataModel and a parameter configuration theta.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.MLE-Tuple{DataModel}","page":"Confidence Regions","title":"InformationGeometry.MLE","text":"MLE(DM::DataModel) -> Vector\n\nReturns the parameter configuration theta_textMLE in mathcalM which is estimated to have the highest likelihood of producing the observed data (under the assumption that the specified model captures the true relationship present in the data). For performance reasons, the maximum likelihood estimate is stored as a part of the DataModel type.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.LogLikeMLE-Tuple{DataModel}","page":"Confidence Regions","title":"InformationGeometry.LogLikeMLE","text":"LogLikeMLE(DM::DataModel) -> Real\n\nReturns the value of the log-likelihood ell when evaluated at the maximum likelihood estimate, i.e. ell(mathrmdata    theta_textMLE). For performance reasons, this value is stored as a part of the DataModel type.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"using InformationGeometry, Plots; gr() # hide\nDS = DataSet([1,2,3.],[4,5,6.5],[0.5,0.45,0.6])\nmodel(x,θ) = θ[1] .* x .+ θ[2]\nDM = DataModel(DS,model)\nMLE(DM), LogLikeMLE(DM)","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"Depending on how the parameters theta enter into the model, the shapes of confidence regions associated with the model may be distorted. For the linearly parametrized model shown above, the 1 sigma and 2 sigma confidence regions form perfect hyperellipses as expected:","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"sols = MultipleConfidenceRegions(DM,1:2)\nVisualizeSols(sols)\n# plot(sols[1],vars=(1,2),label=\"1σ CR\",title=\"Confidence Regions for linearly parametrized model\", xlabel=\"θ[1]\", ylabel=\"θ[2]\") # hide\n# plot!(sols[2],vars=(1,2),label=\"2σ CR\") # hide\n# scatter!([MLE[1]],[MLE[2]],marker=:c,label=\"MLE\") # hide\n# savefig(\"../assets/sols.svg\"); nothing # hide","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"(Image: )","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"For a non-linearly parametrized model, the confidence regions are found to be non-ellipsoidal:","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"model2(x,θ) = θ[1]^3 .* x .+ exp(θ[1] + θ[2])\nDM2 = DataModel(DS,model2)\nsols2 = MultipleConfidenceRegions(DM2,1:2)\nVisualizeSols(sols2)\n#plot(sols2[1],vars=(1,2),label=\"1σ CR\",title=\"Confidence Regions for non-linearly parametrized model\", xlabel=\"θ[1]\", ylabel=\"θ[2]\") # hide\n#plot!(sols2[2],vars=(1,2),label=\"2σ CR\") # hide\n#MLE2 = FindMLE(DM2);  scatter!([MLE2[1]],[MLE2[2]],marker=:c,label=\"MLE\") # hide\n#savefig(\"../assets/sols2.svg\"); nothing # hide","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"(Image: )","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"Specifically in the case of two-dimensional parameter spaces as shown here, the problem of finding the exact boundaries of the confidence regions is turned into a system of ordinary differential equations and subsequently solved using the DifferentialEquations.jl suite. As a result, the boundaries of the confidence regions are obtained in the form of ODESolution objects, which come equipped with elaborate interpolation methods.","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"Since both finding and visualizing exact confidence regions for models depending on more than two parameters (i.e. mathrmdim  mathcalM  2) is more challenging from a technical perspective, the above methods only work for mathrmdim  mathcalM = 2 at this point in time. However, methods which allow for visualizations of confidence regions in arbitrary three-dimensional subspaces of parameter manifolds of any dimension are close to being finished and will follow soon.","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"Various geometric quantities which are intrinsic to the parameter manifold mathcalM can be computed as a result of the Fisher metric g (and subsequent choice of the Levi-Civita connection) such as the Riemann and Ricci tensors and the Ricci scalar R.","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"FisherMetric(::DataModel,::Vector{Float64})\nGeometricDensity(::DataModel,::Vector{Float64})\nChristoffelSymbol(::DataModel,::Vector{Float64})\nRiemann(::DataModel,::Vector{Float64})\nRicci(::DataModel,::Vector{Float64})\nRicciScalar(::DataModel,::Vector{Float64})","category":"page"},{"location":"confidence-regions/#InformationGeometry.FisherMetric-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.FisherMetric","text":"FisherMetric(DM::DataModel, θ::Vector{<:Number})\n\nComputes the Fisher metric g given a DataModel and a parameter configuration theta under the assumption that the likelihood L(mathrmdata    theta) is a multivariate normal distribution.\n\ng_ab(theta) coloneqq -int_mathcalD mathrmd^m y_mathrmdata  L(y_mathrmdata  theta)  fracpartial^2  mathrmln(L)partial theta^a  partial theta^b = -mathbbE bigg( fracpartial^2  mathrmln(L)partial theta^a  partial theta^b bigg)\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.GeometricDensity-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.GeometricDensity","text":"GeometricDensity(DM::DataModel, θ::Vector)\n\nComputes the square root of the determinant of the Fisher metric sqrtmathrmdetbig(g(theta)big) at the point theta.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.ChristoffelSymbol-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.ChristoffelSymbol","text":"ChristoffelSymbol(DM::DataModel, θ::Vector; BigCalc::Bool=false)\nChristoffelSymbol(Metric::Function, θ::Vector; BigCalc::Bool=false)\n\nCalculates the components of the (1,2) Christoffel symbol Gamma at a point theta (i.e. the Christoffel symbol \"of the second kind\") through finite differencing of the Metric. Accurate to ≈ 3e-11. BigCalc=true increases accuracy through BigFloat calculation.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.Riemann-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.Riemann","text":"Riemann(DM::DataModel, θ::Vector; BigCalc::Bool=false)\nRiemann(Metric::Function, θ::Vector; BigCalc::Bool=false)\n\nCalculates the components of the (13) Riemann tensor by finite differencing of the Metric. BigCalc=true increases accuracy through BigFloat calculation.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.Ricci-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.Ricci","text":"Ricci(DM::DataModel, θ::Vector; BigCalc::Bool=false)\nRicci(Metric::Function, θ::Vector; BigCalc::Bool=false)\n\nCalculates the components of the (02) Ricci tensor by finite differencing of the Metric. BigCalc=true increases accuracy through BigFloat calculation.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.RicciScalar-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.RicciScalar","text":"RicciScalar(DM::DataModel, θ::Vector; BigCalc::Bool=false)\nRicciScalar(Metric::Function, θ::Vector; BigCalc::Bool=false)\n\nCalculates the Ricci scalar by finite differencing of the Metric. BigCalc=true increases accuracy through BigFloat calculation.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"Further, studying the geodesics / autoparallels on a manifold can yield enlightening insights about its geometry.","category":"page"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"Score(::DataModel,::Vector{Float64})\nAIC(::DataModel,::Vector{Float64})\nBIC(::DataModel,::Vector{Float64})\nGeodesicDistance(::DataModel,::Vector{Float64},::Vector{Float64})","category":"page"},{"location":"confidence-regions/#InformationGeometry.Score-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.Score","text":"Score(DM::DataModel, θ::Vector{<:Number}; Auto::Bool=false)\n\nCalculates the gradient of the log-likelihood with respect to a set of parameters p. Auto=true uses automatic differentiation.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.AIC-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.AIC","text":"AIC(DM::DataModel, θ::Vector) -> Real\n\nCalculates the Akaike Information Criterion given a parameter configuration theta defined by mathrmAIC = 2  mathrmlength(theta) -2  ell(mathrmdata    theta).\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.BIC-Tuple{DataModel,Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.BIC","text":"BIC(DM::DataModel, θ::Vector) -> Real\n\nCalculates the Bayesian Information Criterion given a parameter configuration theta defined by mathrmAIC = mathrmln(N) cdot mathrmlength(theta) -2  ell(mathrmdata    theta) where N is the number of data points.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/#InformationGeometry.GeodesicDistance-Tuple{DataModel,Array{Float64,1},Array{Float64,1}}","page":"Confidence Regions","title":"InformationGeometry.GeodesicDistance","text":"GeodesicDistance(DM::DataModel,P::Vector{<:Real},Q::Vector{<:Real}; tol::Real=1e-10)\nGeodesicDistance(Metric::Function,P::Vector{<:Real},Q::Vector{<:Real}; tol::Real=1e-10)\n\nComputes the length of a geodesic connecting the points P and Q.\n\n\n\n\n\n","category":"method"},{"location":"confidence-regions/","page":"Confidence Regions","title":"Confidence Regions","text":"To be continued...","category":"page"},{"location":"todo/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"todo/","page":"Contributing","title":"Contributing","text":"If you encounter a bug, feel free to file an issue detailing the problem in contrast to the behaviour you were expecting. Please provide a minimal working example and make sure to specify the particular version of InformationGeometry.jl that was used.\nWhile pull requests are very much welcome, please try to provide detailed docstrings for all non-trivial methods.","category":"page"},{"location":"todo/#Todo:","page":"Contributing","title":"Todo:","text":"","category":"section"},{"location":"todo/","page":"Contributing","title":"Contributing","text":"Allow for non-normal uncertainties in measurements e.g. by interpolating and deriving the Kullback-Leibler divergence over a domain\nAllow for uncertainties in the conditions (i.e. x-values) by rewriting the DataSet type\nParallelism: Improve support for parallel computations of geodesics, curvature tensors and so on\nEmploy importance sampling for Monte Carlo computations\nImprove visualization capabilities for high-dimensional models\nImprove documentation, docstrings and provide more examples\nStandardize the user-facing keyword arguments\nProvide performance benchmarks for InformationGeometry.jl\nUse IntervalArithmetic.jl and IntervalOptimisation.jl for rigorous guarantees on inference results?","category":"page"},{"location":"datamodels/#Providing-Datasets","page":"Providing Data and Models","title":"Providing Datasets","text":"","category":"section"},{"location":"datamodels/","page":"Providing Data and Models","title":"Providing Data and Models","text":"Typically, one of the most difficult parts of any data science problem is to bring the data into a form which lends itself to the subsequent analysis. This section aims to describe the containers used by InformationGeometry.jl to store datasets and models in detail.","category":"page"},{"location":"datamodels/","page":"Providing Data and Models","title":"Providing Data and Models","text":"The data itself is stored using the DataSet container.","category":"page"},{"location":"datamodels/","page":"Providing Data and Models","title":"Providing Data and Models","text":"DataSet","category":"page"},{"location":"datamodels/#InformationGeometry.DataSet","page":"Providing Data and Models","title":"InformationGeometry.DataSet","text":"The DataSet type is a versatile container for storing data. Typically, it is constructed by passing it three vectors x, y, sigma where the components of sigma quantify the standard deviation associated with each y-value. Alternatively, a full covariance matrix can be supplied for the ydata instead of a vector of standard deviations. The contents of a DataSet DS can later be accessed via xdata(DS), ydata(DS), sigma(DS).\n\nExamples:\n\nIn the simplest case, where all data points are mutually independent and have a single x-component and a single y-component each, a DataSet consisting of four points can be constructed via\n\nDataSet([1,2,3,4],[4,5,6.5,7.8],[0.5,0.45,0.6,0.8])\n\nor alternatively by\n\nDataSet([1,2,3,4],[4,5,6.5,7.8],Diagonal([0.5,0.45,0.6,0.8].^2))\n\nwhere the diagonal covariance matrix in the second line is equivalent to the vector of uncertainties supplied in the first line.\n\nMore generally, if a dataset consists of N points where each x-value has n many components and each y-value has m many components, this can be specified to the DataSet constructor via a tuple (Nnm) in addition to the vectors x, y and the covariance matrix. For example:\n\nX = [0.9, 1.0, 1.1, 1.9, 2.0, 2.1, 2.9, 3.0, 3.1, 3.9, 4.0, 4.1]\nY = [1.0, 5.0, 4.0, 8.0, 9.0, 13.0, 16.0, 20.0]\nCov = Diagonal([2.0, 4.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0])\ndims = Tuple([4,3,2])\nDS = DataSet(X,Y,Cov,dims)\n\n\n\n\n\n","category":"type"},{"location":"datamodels/","page":"Providing Data and Models","title":"Providing Data and Models","text":"To complete the specification of an inference problem, a model function which is assumed to model the relationship which is inherent in the takes an x-value and a parameter configuration theta must be added.","category":"page"},{"location":"datamodels/","page":"Providing Data and Models","title":"Providing Data and Models","text":"DataModel","category":"page"},{"location":"datamodels/#InformationGeometry.DataModel","page":"Providing Data and Models","title":"InformationGeometry.DataModel","text":"In addition to storing a DataSet, a DataModel also contains a function model(x,θ) and its derivative dmodel(x,θ) where x denotes the x-value of the data and θ is a vector of parameters on which the model depends. Crucially, dmodel contains the derivatives of the model with respect to the parameters θ, not the x-values. For example\n\nDataSet([1,2,3,4],[4,5,6.5,7.8],[0.5,0.45,0.6,0.8])\nmodel(x,θ::Vector) = θ[1] .* x .+ θ[2]\nDM = DataModel(DS,model)\n\nIn order for all methods to perform as expected, the output of the model must always be in the form of a vector of numbers (if the y-values have more than one component or the model is evaluated on more than one point).\n\nIf a DataModel is constructed as shown above, the gradient of the model with respect to the parameters θ (i.e. its \"Jacobian\") will be calculated using automatic differentiation. Alternatively, an explicit analytic expression for the Jacobian can be specified by hand:\n\nfunction dmodel(x,θ::Vector)\n   J = Array{Float64}(undef, length(x), length(θ))\n   @. J[:,1] = x        # ∂(model)/∂θ₁\n   @. J[:,2] = 1.       # ∂(model)/∂θ₂\n   return J\nend\nDM = DataModel(DS,model,dmodel)\n\nThe output of the Jacobian must be a matrix whose columns correspond to the partial derivatives with respect to different components of θ and whose rows correspond to evaluations at different values of x.\n\nThe DataSet contained in a DataModel DM can be accessed via DM.Data, whereas the model and its Jacobian can be used via DM.model and DM.dmodel respectively.\n\n\n\n\n\n","category":"type"},{"location":"basics/#Introduction-to-Information-Geometry","page":"Basics of Information Geometry","title":"Introduction to Information Geometry","text":"","category":"section"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"For a more detailed discussion of information geometry, see e.g. my Master's Thesis, this paper or this book.","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"Essentially, information geometry is a combination of the mathematical disciplines of differential geometry and probability theory. The main idea is to rephrase statistical problems in such a way that they can be given a geometric interpretation.","category":"page"},{"location":"basics/#Information-Divergences","page":"Basics of Information Geometry","title":"Information Divergences","text":"","category":"section"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"In information theory, the dissimilarity between two probability distributions p(x) and q(x) is generally quantified using so-called information divergences, which are positive-definite functionals. The most popular choice of information divergence is given by the Kullback-Leibler divergence D_textKLpq defined by","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"D_textKLpq coloneqq int mathrmd^m y  p(y)  mathrmln bigg( fracp(y)q(y) bigg) = mathbbE_p biggmathrmlnbigg( fracpq bigg) bigg","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"Intuitively, the Kullback-Leibler divergence corresponds to the relative increase in Shannon entropy (i.e. loss of information) that is incurred by approximating the distribution p(x) through q(x). In addition to its tangible information-theoretic interpretation, the Kullback-Leibler divergence has the following desirable properties:","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"reparametrization invariance with respect to the random variable over which the distributions are integrated,\napplicable between any two probability distributions with common support, e.g. a chi^2-distribution and a Poisson distribution or a normal and a Cauchy distribution.","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"On the other hand, the disadvantages of using information divergences such as the Kullback-Leibler divergence to measure the dissimilarity of distributions are:","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"they are typically not symmetric, i.e. D_textKLpq neq D_textKLqp\nthey usually do not satisfy a triangle inequality","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"wherefore they do not constitute distance functions (i.e. metric functions) on the underlying space of probability distributions.","category":"page"},{"location":"basics/#The-Fisher-Metric","page":"Basics of Information Geometry","title":"The Fisher Metric","text":"","category":"section"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"In practical applications, one is often particularly interested in spaces of probability distributions which form a single overarching family and can be parametrized using a parameter configuration theta in mathcalM where mathcalM constitutes a smooth manifold. Accordingly, any two members p(ytheta_1) and p(ytheta_2) of this family can be compared using e.g. the Kullback-Leibler divergence D_textKLbigp(theta_1)p(theta_2)big via the formula given above.","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"While the Kullback-Leibler divergence D_textKLpq does not constitute a proper distance function on mathcalM, it can be expanded in Taylor series around theta_textMLE in terms of its derivatives. The zeroth order of this expansion vanishes due to the definiteness of the Kullback-Leibler divergence (i.e. D_textKLqq = 0 for all distributions q). Similarly, the first order vanishes since the expectation of the components of the score is nil. Thus, the second order approximation of the Kullback-Leibler divergence is completely determined by its Hessian, which can be computed as","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"g_ab(theta) coloneqq biggfracpartial^2partial psi^a  partial psi^b  D_textKL bigp(ytheta)  p(ypsi) big bigg_psi = theta\n=  = -mathbbE_pbiggfracpartial^2  mathrmln(p)partial theta^a  partial theta^bbigg =  = mathbbE_pbiggfracpartial  mathrmln(p)partial theta^a fracpartial  mathrmln(p)partial theta^bbigg","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"where it was assumed that the order of the derivative operator and the integration involved in the expectation value can be interchanged.","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"The Hessian of the Kullback-Leibler divergence is typically referred to as the Fisher information matrix. Moreover, since it can be shown that the Fisher information is not only positive-definite but also exhibits the transformation behaviour associated with a (02)-tensor field, it can therefore be used as a Riemannian metric on the parameter manifold mathcalM.","category":"page"},{"location":"basics/","page":"Basics of Information Geometry","title":"Basics of Information Geometry","text":"Clearly, the Riemannian geometry induced on mathcalM by the Fisher metric is ill-equipped to faithfully capture the behaviour of the Kullback-Leibler divergence in its entirety (e.g. its asymmetry). Nevertheless, this Riemannian approximation already encodes many of the key aspects of the Kullback-Leibler divergence and additionally benefits from the versatility and maturity of the differential-geometric formalism. Therefore, the Fisher metric offers a convenient and powerful tool which can be used to study statistical problems in a coordinate invariant setting which focuses on intrinsic properties of the parameter manifold.","category":"page"},{"location":"#InformationGeometry","page":"Getting Started","title":"InformationGeometry","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"This is the documentation of InformationGeometry.jl, a Julia package for differential-geometric analyses of statistical problems.","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Build Status\n(Image: travis) (Image: appveyor) (Image: codecov)","category":"page"},{"location":"#Main-Uses","page":"Getting Started","title":"Main Uses","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"maximum likelihood estimation\nconstruction and visualization of exact confidence regions\nefficient calculation of Kullback-Leibler divergences\ncomputation of geometric quantities such as geodesics and curvature on the parameter manifold","category":"page"},{"location":"#Installation","page":"Getting Started","title":"Installation","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"As with any Julia package, InformationGeometry.jl can be added from the Julia terminal via","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"julia> ] add InformationGeometry","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"or alternatively by","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"julia> using Pkg; Pkg.add(\"InformationGeometry\")","category":"page"},{"location":"kullback-leibler/#Kullback-Leibler-Divergences","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"","category":"section"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"Using the Distributions type provided by Distributions.jl, the KullbackLeibler method offers a convenient way of computing the Kullback-Leibler divergence between distributions. In several cases an analytical expression for the Kullback-Leibler divergence is known. These include: (univariate and multivariate) Normal, Cauchy, Exponential, Weibull and Gamma distributions.","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"Furthermore, for distributions over a one-dimensional domain where no analytic result is known, KullbackLeibler rephrases the integral in terms of an ODE and employs an efficient integration scheme from the DifferentialEquations.jl suite. For multivariate distributions, Monte Carlo integration is used.","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"Examples of use:","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"KullbackLeibler(Cauchy(1.,2.4),Normal(-4,0.5),HyperCube([-100,100]);Carlo=false,tol=1e-12)\nKullbackLeibler(MvNormal([0,2.5],diagm([1,4.])),MvTDist(1,[3,2],diagm([2.,3.])),HyperCube([[-50,50],[-50,50]]); N=Int(1e8))","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"In addition, it is of course also possible to input generic functions, whose positivity and normalization should be ensured by the user.","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"KullbackLeibler(::Function,::Function)","category":"page"},{"location":"kullback-leibler/#InformationGeometry.KullbackLeibler-Tuple{Function,Function}","page":"Kullback-Leibler Divergences","title":"InformationGeometry.KullbackLeibler","text":"KullbackLeibler(p::Function,q::Function,Domain::HyperCube=HyperCube([[-15,15]]); tol=2e-15, N::Int=Int(3e7), Carlo::Bool=(Domain.dim!=1))\n\nComputes the Kullback-Leibler divergence between two probability distributions p and q over the Domain. If Carlo=true, this is done using a Monte Carlo Simulation with N samples. If the Domain is one-dimensional, the calculation is performed without Monte Carlo to a tolerance of ≈ tol.\n\nD_textKLpq coloneqq int mathrmd^m y  p(y)  mathrmln bigg( fracp(y)q(y) bigg)\n\n\n\n\n\n","category":"method"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"For example, the Kullback-Leibler divergence between a Cauchy distribution with mu=1 and s=2 and a normal (i.e. Gaussian) distribution with mu=-4 and sigma=12 can be calculated via:","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"using InformationGeometry # hide\nusing LinearAlgebra, Distributions\nKullbackLeibler(Cauchy(1.,2.),Normal(-4.,0.5),HyperCube([-100,100]); Carlo=false,tol=1e-12)","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"Specifically, the keyword arguments used here numerically compute the divergence over the domain -100100 to an accuracy of 10^-12.","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"The domain of the integral involved in the computation of the divergence is specified using the HyperCube datatype, which stores a cuboid region in N dimensions as a vector of intervals.","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"HyperCube","category":"page"},{"location":"kullback-leibler/#InformationGeometry.HyperCube","page":"Kullback-Leibler Divergences","title":"InformationGeometry.HyperCube","text":"The HyperCube type has the fields vals::Vector{Vector}, which stores the intervals which define the hypercube and dim::Int, which gives the dimension. Overall it just offers a convenient and standardized way of passing domains for integration or plotting between functions without having to check that these domains are sensible every time. Examples for constructing HyperCubes:\n\nHyperCube([[1,3],[pi,2pi],[-500.0,100.0]])\nHyperCube([[-1,1]])\nHyperCube([-1,1])\nHyperCube(LowerUpper([-1,-5],[0,-4]))\nHyperCube(collect([-7,7.] for i in 1:3))\n\nThe HyperCube type is closely related to the LowerUpper type and they can be easily converted into each other. Examples of quantities that can be computed from and operations involving a HyperCube object X:\n\nCubeVol(X)\nTranslateCube(X,v::Vector)\nCubeWidths(X)\n\n\n\n\n\n","category":"type"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"Furthermore, the Kullback-Leibler divergence between multivariate distributions can be computed for example by","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"KullbackLeibler(MvNormal([0,2.5],diagm([1,4.])),MvTDist(1,[3,2],diagm([2.,3.])),HyperCube([[-50,50],[-50,50]]); N=Int(5e6))","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"where it now becomes necessary to employ Monte Carlo schemes. Specifically, the keyword argument N now determines the number of points where the integrand is evaluated over the domain -5050 times -5050.","category":"page"},{"location":"kullback-leibler/","page":"Kullback-Leibler Divergences","title":"Kullback-Leibler Divergences","text":"So far, importance sampling has not been implemented for the Monte Carlo integration. Instead, the domain is sampled uniformly.","category":"page"}]
}

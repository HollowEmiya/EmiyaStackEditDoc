# New 大气散射

## Transmittance 透光率

光线在大气层中从一点 **$\mathrm{p}$** 到达一点 **$\mathrm{q}$** ，由于空气中粒子的存在可能会发生光的能量被吸收或者散射到别的地方去。  
因此到达 **$\mathrm{q}$** 的光线只是 **$\mathrm{p}$** 的一部分，这部分光存留多少取决于光的波长，被称为 [Transmittance]([Transmittance - Wikipedia](https://en.wikipedia.org/wiki/Transmittance)) (透光率)。

人话：光经过一段距离还剩多少。  
Transmittance 描述的是这一段和方向无关的距离，因为这段大气在这里，你正着进来和反着进来都会衰减一样的。

### 计算

定义物理单位，这部可以不做，做了可以方便修改所有计算的物理尺度(我是懒狗所以我不做)

~~~C++
const float m = 1.0f;
const float km = 1000.0f * m;
const float m2 = m * m;

const float rad = 1.0f;
~~~

假设有共线三点顺序为 A->B->C，且都在大气层中, 按这种顺序 AC 之间的 Transmittance 是 AB 的 Transmittance 和 BC 的 Transmittance 的乘积。
$$
T(\overline{AC})=T(\overline{AB})*T(\overline{BC})
$$
 所以，下图 pq 之间的 transmittance 是，有射线[p,q) 和大气顶部(大气层)或大气底部(大气层底部这里在地球上可以认为是地球表面)交点 i，用 $\overline{pi}$ 的 transmittance 除 $\overline{qi}$ 的 transmittance, 即 $T(\overline{pq})=\frac{T(\overline{pi})}{T(\overline{qi})}$.(当然如果[p,q]和地面相交那 Transmittance 就是 0)

*<del>什么顶部底部，tmd不说人话。</del>*

这里只模拟大气散射，所以对一个星球我们只关心 大气 和 非大气部分，大气底部的球体大小用 atmosphere_bottom_radius 定义，而 大气顶部即 非大气球体+大气厚度 用 atmosphere_top_radius 定义。
下图的最下面的曲线表示的就是非大气部分——atmosphere_bottom_radius, 最外边的是 atmosphere_top_radius.

![enter image description here](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/computation-0.png?raw=true)

还有 pq 的 transmittance 和 qp 的 transmittance 是一样的，这和方向无关。  
所以计算任意两点间的 Transmittance 我们只需要知道 一点 **p** 和与大气顶部的交点 **i**。  
(这也是为什么可以用预计算进行大气模拟，因为我们可以对 Transmittance 做无限高的预计算，反正不会影响我们的实时渲染)
而 Transmittance 只受两个变量的影响，一个是 上图 p 到地心 的距离 $$r=\vert\vert op\vert\vert$$ 和 光线方向 和天顶的夹角，就是 射线 和 $\vec{up}$方向的夹角。
$\cos(\theta)=\frac{\vec {op}\cdot\vec{pi}}{\vert\vec{op}\vert\vert\vec{pi}\vert}$, 所以我们得先知道 $\vec{pi}$ 的长度，也就是 p 到大气外部的交点，并且要知道什么时候，线段 $\vec{pi}$ 和地面相交。

这样就能把 Transmittance 映射到 2维纹理上。

假设一点在沿 $\vec{pi}$ 方向上，且距离 p 点距离为 d，则有坐标(这里认为行星中心是0,0) $[d\sqrt{1-\mu^2},r+d\mu]$.  ，所以关于点 i ，我们能得到 $\vert\vert pi\vert\vert^2+2r\mu\vert\vert pi\vert\vert+r^2=r^2_{top}$
先计算第一种情况，就是交点在大气部分

~~~C++
// planetAtmoRadius is the planet radius + Atmosphere radius
float DistanceToTopAtmosphereBoundary(float planetAtmoRadius, float r, float cosTheta)
{
	assert(r <= planetAtmoRadius)
	assert(-1.0 <= cosTheta && cosTheta <= 1.0)
	// b*b - 4ac
	float disSquare = r*r *(cosTheta*cosTheta - 1.0) + 
					planetAtmoRadius * planetAtmoRadius;
	return -r*cosTheta + sqrt(disSquare);
}
~~~

为什么是+，首先 p 点在 大气球层的内部，而 i 是 沿 $\theta$ 方向前进距离 d 的一点，所以 d 一定要保证大于零，故而是 +，如果是 - 就是反方向了。

现在求 p 和大气底部的交点，(这里的计算是假设一定有交点)

~~~c++
float DistanceToBottomAtmosphereBoundary(float atmoBottomRadius, float r, float cosTheta)
{
	assert(r >= atmoBottomRadius)
	assert(-1.0 <= cosTheta && cosTheta <= 1.0)
	// b*b - 4ac
	float disSquare = r*r *(cosTheta*cosTheta - 1.0) + 
					atmoBottomRadius * atmoBottomRadius;
	return -r*cosTheta - sqrt(disSquare);
}
~~~

这里为什么又是 - ，
因为假设一定相交，所以要取距离近的交点，所以是 - 。

让 [**p**,**i**] 和非大气球有交点，要保证 $d^2+2r\mu d+r^2=r^2_{bottom}$ 有非负解。  
$r^2(\mu^2-1)+r^2_{bottom} \ge0$ 函数可以写作

~~~C++
bool RayIntersectsGround(float planetRadius, float r, float cosTheta)
{
	return cosTheta < 0.0 && r*r + (cosTheta*cosTheta - 1) +
		planetRadius*planetRadius >= 0.0f;
}
~~~

现在计算 [**p**,**i**] 之间的 transmittance，根据定义和 [比尔-朗伯定理](https://en.wikipedia.org/wiki/Beer–Lambert_law), 这受影响于 [**p**,**i**] 之间的空气分子浓度的积分和 气溶胶的数量密度积分 以及 吸收光的空气粒子(臭氧)数量浓度积分。  
这三个积分形式一致，当 [**p**,**i**] 不和非大气球相交时，可以使用[梯形公式](https://en.wikipedia.org/wiki/Trapezoidal_rule)计算  
(这里计算是认为  [**p**,**i**] 和大气球体相交)  

~~~C++
// 先看一下描述密度分布的模型吧
// width 就是该 气体层级 的起始高度
// 'exp_term' * exp('exp_scale' * h) + 'linear_term' * h + 'constant_term'
// h 是距离非大气球壳的高度
struct DensityProfileLayer
{
	// 
	DensityProfileLayer():
		DensityProfileLayer(0.0, 0.0, 0.0, 0.0, 0.0){}
	
	DensityProfileLayer(float width, float exp_term, float exp_scale, 
		float linear_term, float constant_term)
		: width(width), exp_term(exp_term), exp_scale(exp_scale),
        linear_term(linear_term), constant_term(constant_term) 
    {
  	}

	float width;
	float exp_term;
	float exp_scale;
	float linear_term;
	float constant_term;
};

struct DensityProfile
{
	DensityProfileLayer layers[2];
};

// 这里拿 臭氧层举例吧……
// 臭氧层是中间厚两边薄的。
// 这里 DensityProfile 分两个 layer 也是为了这个
// layer 0 表示 [0,25] 的部分
// layer 1 表示 25 km 以上部分
//layers[0] 0.0 * exp(0.0 * h) + 1.0/15.0 * h - 2.0/3.0
//layers[1] 0.0 * exp(0.0 * h) - 1.0/15.0 * h + 8.0/3.0
DensityProfile ozoneProfile.layers[0] = DensityProfileLayer(
        25.0 * km, 0.0, 0.0 / km, 1.0 / (15.0 * km), -2.0 / 3.0);
    ozoneProfile.layers[1] = DensityProfileLayer(
        0.0 * km, 0.0, 0.0 / km, -1.0 / (15.0 * km), 8.0 / 3.0);

float GetLayerDensity(DensityProfileLayer layer, float altitude)
{
	float density = layer.exp_term * exp(layer.exp_scale * altitude) + 
		layer.linear_term * altitude + layer.constant_term;
    return density;
}

float GetProfileDensity(DensityProfile profile, float altitude)
{
	// 以臭氧为例，低于 25 km，采用layer0，否则采用layer2
    // 也不一定是臭氧，可能是所有 abort 粒子……不知道
	return altitude < profile.layer[0].width ?
		GetLayerDensity(profile.layers[0], altitude) :
		GetLayerDensity(profile.layers[1], altitude);
}

float ComputeOpticalLengthToAtmosphereBoundary(
	AtmosphereParameters atmosphere, DensityProfile profile, float r, float cosTheta
)
{
	assert(r <= atmosphere.topRadius && r >= atmosphere.bottomRadius)
	assert(-1.0 <= cosTheta && cosTheta <= 1.0)
	
	int SAMPLE_COUNT = 500;
	float dx = DistanceToTopAtmosphereBoundary(atmosphere.topRadius, r, cosTheta) / SAMPLE_COUNT;
	float res = 0;
	for(int i = 0; i <= SAMPLE_COUNT; ++i)
	{
		float d_i = i * dx;
		// (r + d_i*cosTheta)^2 + (d_i*sinTheta)^2 = (r_i)^2
		float r_i = sqrt(r*r + 2.0*r*d_i*cosTheta + d_i*d_i);
		float y_i = GetProfileDensity(profile, r_i - atmosphere.bottomRadius);
		float weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
		res += y_i * weight_i * dx;
	}
	return result;
}
~~~
[Nasa Ozone Watch: Ozone facts](https://ozonewatch.gsfc.nasa.gov/facts/SH.html)
臭氧层是中间厚两边薄的。

> The total mass of ozone in the atmosphere is about 3 billion metric tons. That may seem like a lot, but it is only 0.00006 percent of the atmosphere. The peak concentration of ozone occurs at an altitude of roughly 32 kilometers (20 miles) above the surface of the Earth. At that altitude, ozone concentration can be as high as 15 parts per million (0.0015 percent).
> ![enter image description here](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/ozone_concentration_graph.gif?raw=true)


这样就能计算 p，i 之间的 transmittance 了，这里假设pi段不会和地面相交
~~~c++
struct AtmoSphereParameters
{
    ...
    DensityProfile rayleigh_density;
    DensityProfile mie_density;
    DensityProfile absorption_density;
	
    float3 rayleigh_scattering;
    float3 mie_scattering;
    float3 absorption_extinction;
    
    float bottomRadius;	// planetRaiuds
    float topRadius;	// planetAtmoRadius
    ...

	// The cosine of the maximum Sun zenith angle 
	// for which atmospheric scattering
	// must be precomputed 
	// (for maximum precision, use the smallest Sun zenith
	// angle yielding negligible sky light radiance values.
	// For instance, for the Earth case, 102 degrees 
	// is a good choice - yielding mu_s_min = -0.2).

	// 余弦值的最大太阳天顶角，需要对大气散射进行预计算
	//（为了最大精度，请使用使得天空光辐射值可以忽略的最小太阳天顶角。
	// 例如，对于地球的情况，选择102度是一个不错的选择，
	// 对应的μ_s_min为-0.2）。
	float mu_s_min;
};

float3 ComputeTransmittanceToTopAtmosphereBoundary(
	AtmoSphereParameters atmosphere, float r, float cosTheta
)
{
    assert( r >= atmosphere.bottomRaiud && r <= planetAtmoRadius);
    assert( cosTheta <= 1 && cosTheta >= -1);
    
    return exp( atmosphere.rayleight_scattering *
              ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere,
				atmosphere.rayleigh_density, r, cosTheta) + 
           atmosphere.mie_scattering * 
           	ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere,
				atmosphere.mie_density, r, cosTheta) + 
         	atmosphere.absorption_density *
               ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere,
				atmosphere.absorption_extinction, r, cosTheta) + 
               );
}
~~~

### Precomputation

上述函数的计算成本相当高，而且需要大量的计算来进行单次和多次散射的计算。幸运的是，这个函数只依赖于两个参数，并且是相当平滑的，因此我们可以通过在一个小的2D纹理中预计算它来优化其计算。

为此，我们需要一个函数参数（r, μ）和纹理坐标（u, v）之间的映射，以及其反向映射，因为这些参数的单位和值范围不同。即使是在这种情况下，由于纹理采样在纹理元素中心，将函数 f 从 [0,1] 区间存储在大小为 n 的纹理中会在 0.5/n、1.5/n、... (n−0.5)/n 处进行采样。因此，这个纹理只会在域边界（0 和 1）处给我们外推的函数值。为了避免这种情况，我们需要在纹理元素0的中心存储 f(0) ，并在纹理元素 n−1 的中心存储 f(1) 。这可以通过以下从值 x 在 [0,1] 区间到纹理坐标 u 在 [0.5/n,1−0.5/n] 区间的映射及其反向映射来实现：


> For this we need a mapping between the function parameters (r,μ) and the texture coordinates (u,v), and vice-versa, because these parameters do not have the same units and range of values. And even if it was the case, storing a function f from the [0,1] interval in a texture of size n would sample the function at 0.5/n, 1.5/n, ... (n−0.5)/n, because texture samples are at the center of texels. Therefore, this texture would only give us extrapolated function values at the domain boundaries (00 and 11). To avoid this we need to store f(0) at the center of texel 0 and f(1) at the center of texel n−1. This can be done with the following mapping from values x in [0,1] to texture coordinates u in [0.5/n,1−0.5/n] - and its inverse:

~~~c++
// 从 [0,1] 到防止插值的空间，进行采样
float GetTextureCoordFromUnitRange(float x, int tex_size)
{
    return 0.5/ tex_size + x * (1.0 - 1.0/tex_size);
}

// 从防插值空间变换到 [0,1] 进行数学计算
float GetUnitRangeFromTextureCoord(float u, int tex_size)
{
    return (u - 0.5 / tex_size) / (1.0 - 1.0 / tex_size);
}
~~~

沿用这个思路我们可以定义从 $(r,\mu)$ 到 $(u,v)$ 的映射，这样就能避免 纹理采样时在边界值0或1时采样错误。在 [2008 Eric 的论文 Precomputed Atmospheric Scattering](https://inria.hal.science/inria-00288758/en) 中这种映射关系使用了很多只适用于地球大气的参数。这里将使用一种更通用的映射方式，适用于任何大气，不过仍然需要在视界线(地平线、 horizon)依然要提升采样率。这里的映射方式基于 使用 4D texture 的 [paper](https://inria.hal.science/inria-00288758/en) 对 r 使用相同的映射，对 cosTheta 进行一些改动 (只考虑视线不和地面相交的情况)。具体来讲将 cosTheta 映射到一个 $x_\mu$，$x_\mu$ 受到 距离大气球层距离 d 的影响，且取值范围[0,1]，
$$
d_{min}=r_{top}-r\\
d_{max}=\rho+H\\
\rho=\sqrt {r^2-bottom\_radius^2}\\
H=\sqrt{top\_radius^2-bottom\_radius^2}
$$

![enter image description here](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/precomputation-0.png?raw=true)

~~~C++
float2 GetTransmittanceTextureUvFrom_RAndCosTheta(
    AtmosphereParameters atmosphere, float r, float cosTheta)
{
    assert( r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
    assert( cosTheta >= -1.0 && cosTheta <= 1.0);
    
    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
		atmosphere.bottom_radius * atmosphere.bottom_radius);
	
	float rho =
		sqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float d = DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    return float2(
        GetTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_TEXTURE_WIDTH),
        GetTextureCoordFromUnitRange(x_r, TRANSMITTANCE_TEXTURE_HEIGHT)
    );
}
~~~

![enter image description here](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/atmosphere.png?raw=true)
从这里的 d_min 和 d_max 可以看出这里所有的计算值都是在大气层中且不会和地面相交  

逆 uv 变换为：

~~~C++
void GetRMuFromTransmittanceTextureUv(AtmosphereParameters atmosphere,
    float2 uv, out float r, out float cosTheta) {
	assert(uv.x >= 0.0 && uv.x <= 1.0);
    assert(uv.y >= 0.0 && uv.y <= 1.0);
    float x_cosTheta = GetUnitRangeFromTextureCoord(uv.x, TRANSMITTANCE_TEXTURE_WIDTH);
    float x_r = GetUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT);

    float  H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
      atmosphere.bottom_radius * atmosphere.bottom_radius);

    float rho = H * x_r;
    r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);

    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float d = d_min + x_cosTheta * (d_max - d_min);
    cosTheta = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
    mu = Clamp(mu,-1.0,1.0);
}
~~~

![enter image description here](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/atmosphere_%CE%B8.png?raw=true)
[atmosphere - GeoGebra](https://www.geogebra.org/geometry/ngykf3sz)

现在就能够用片元着色器计算 transmittance 了
~~~c++
float3 ComputeTransmittanceToTopAtmosphereBoundaryTexture(
	AtmosphereParameters atmosphere, float2 frag_coord
)
{
    const float2 TRANSMITTANCE_TEXTURE_SIZE =
      vec2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
  float r;
  float cosTheta;
  GetRMuFromTransmittanceTextureUv(
      atmosphere, frag_coord / TRANSMITTANCE_TEXTURE_SIZE, r, cosTheta);
  return ComputeTransmittanceToTopAtmosphereBoundary(atmosphere, r, cosTheta);
}
~~~

### LookUp

有了上述方法所做的 预计算tex，我们就能得到大气球体中任意一点到大气层顶部之间的 transmittance(这里假设前提是这两点不会和地面相交)

> 其实已经可以算任意两点且此两点间的了……

~~~c++
float3 GetTransmittanceToTopAtmosphereBoudary(
	AtmosphereParameters atmosphere,
    TrasmittanceTexture transmittance_texture,
    float r, float cosTheta
)
{
    assert( r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
    float2 uv = GetTransmittanceTextureUvFrom_RAndCosTheta(atmosphere, r, cosTheta);
    return DimensionlessSpectrum(texture(transmittance_texture, uv));
}
~~~

c++设有一点p在大气中，沿和天顶夹角 θ方向和大气顶部交于 i，该路径上有一点 q，q 距离 p 为 d。
则 $r_d=||oq||=\sqrt{d^2+2dr\cos\theta+r^2}$, $\cos\theta_d={\vec{oq}\cdot\vec{pi}}\,/\,{||\vec{oq}||*||\vec{pi}||}=(r\cos\theta+d)/r_d$ 
有了 q 点的 r 和 μ，我们用之前得到的 transmittance_tex 就能够计算大气中 p, q 两点之间的 transmittance 了。  
![enter image description here](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/atmosphere_thetad.png?raw=true)
[atmosphere_rd - GeoGebra](https://www.geogebra.org/geometry/z8a3uxa3)
要计算pq的 transmittance 只需要用 $T(\vec{pq})=T(\vec{pi}) / T(\vec{qi})$

~~~c++
float3 GetTransmittance(
	AtmosphereParameters atmosphere, Texture transmittance_texture,
	float r, float cosTheta, float d, bool ray_r_cosTheta_intersects_ground
)
{
    assert( r >= atmosphere.bottom_radius && r <= 	atmosphere.top_radius);
    assert(cosTheta >= -1.0 && cosTheta <= 1.0);
    assert(d >= 0.0 * m);

    float r_d = ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * cosTheta * d + r * r));
    float mu_d = ClampCosine((r * cosTheta + d) / r_d);

    if (ray_r_cosTheta_intersects_ground) {
        // cross ground T(qp)=T(qi')/T(pi')
    return min(
        GetTransmittanceToTopAtmosphereBoundary(
            atmosphere, transmittance_texture, r_d, -cosTheta_d) /
        GetTransmittanceToTopAtmosphereBoundary(
            atmosphere, transmittance_texture, r, -cosTheta),
        float3(1.0,1.0,1.0));
    } else {
        // p,q,i all in atmosphere, dont cross ground
    return min(
        GetTransmittanceToTopAtmosphereBoundary(
            atmosphere, transmittance_texture, r, cosTheta) /
        GetTransmittanceToTopAtmosphereBoundary(
            atmosphere, transmittance_texture, r_d, cosTheta_d),
        float3(1.0,1.0,1.0));
    }
}
~~~

c++如果以 $r$ 和 $\mu即\cos\theta$ 来定义射线和地球相交，这里的 `ray_r_cosTheta_intersects_ground` 应该为`true`。这里不会用 `RayIntersectsGround` 计算是否相交，因为这样计算的话，当射线非常接近地平线(视界线)时由于浮点数的精度问题会导致计算结果错误。[^0]后面我们有更好的方法来计算射线和地球相交。  
为什么和地面相交要反着算：
![Transmittance](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/Transmittance_pq.png?raw=true)  
[Transmittance_Pq](https://www.geogebra.org/geometry/tuahqrt2)  
最后我们需要大气中某一点到太阳的 Transmittance。太阳不是点光源，所以这是对太阳球盘的积分。  
我们可以认为在太阳球盘的积分是一个常数，除非在地平线(视界线)以下为零的这种情况。  
所以，对太阳的 Transmittance 可以由 `GetTransmittanceToTopAtmosphereBoundary` 计算，然后乘太阳球盘在地平线上的比例。  
当太阳从地面升起，即太阳天顶角 $\theta_s$ 大于 地平线天顶角 $\theta_h$ 加上 太阳角半径 $\alpha_s$，即 $\theta_s>\theta_h+\alpha_s$ ,这个比例从 0 到 1。等价于 $[\mu_s=\cos\theta_s]>[\;\cos(\theta_h+\alpha_s)\approx\;\cos\theta_s-\alpha_s\sin\theta_h]$。  
而当太阳落山时，这个比例从 1 到 0，即  即太阳天顶角 $\theta_s$ 小于 地平线天顶角 $\theta_h$ 减去 太阳角半径 $\alpha_s$，即 $\theta_s<\theta_h-\alpha_s$ 。等价于 $[\mu_s=\cos\theta_s]<[\;\cos(\theta_h-\alpha_s)\approx\;\cos\theta_s+\alpha_s\sin\theta_h]$。  
在这两者之间，可见的太阳圆盘部分的变化近似于一个平滑步长(可以通过绘制圆形[弓形](https://en.wikipedia.org/wiki/Circular_segment)的面积作为其[矢状图](https://en.wikipedia.org/wiki/Sagitta_(geometry))的函数来验证这一点)。因此，由于 $sin\theta_h=r_{bottom}/r$ (不明白就自己画一个图，做p和planet的切线，再过圆心O做该切线的垂线)，我们可以用以下函数近似太阳的透射率:

~~~C++
float3 GetTransmittanceToSun(
	AtmosphereParameters atmosphere,
	Texture transmittance_texture,
    float r, float cosTheta_s
)
{
    float sin_theta_h = atmosphere.bottom_radius / r;
    float cos_theta_h = -sqrt(max(1.0 - sin_theta_h * sin_theta_h, 0.0));
    return GetTransmittanceToTopAtmosphereBoundary(
          atmosphere, transmittance_texture, r, cosTheta_s) *
      smoothstep(-sin_theta_h * atmosphere.sun_angular_radius / rad,
                 sin_theta_h * atmosphere.sun_angular_radius / rad,
                 cosTheta_s - cos_theta_h);
}
~~~

![](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/Sun.png?raw=true)  
$$
\begin{aligned}
In\;fact\;\alpha_s\;is\;due\;to\;r,but\;the\;dis\;from\;earth\;to\;sun\;is\;very\;far.\\
Here\; we\;make\;\alpha_s\;constant.\;Whatever\;this\;is\;not\;the\;point.\hspace{35px}\\
\theta_s <\theta_h-\alpha_s\hspace{40px}The\; rate\;is\;1.\hspace{275px}\\
\theta_s >\theta_h+\alpha_s\hspace{40px}The\; rate\;is\;0.\hspace{275px}\\
\theta_h+\alpha_s>\theta_s >\theta_h-\alpha_s\hspace{40px}The\; rate\;is\;[0,1].\hspace{158px}\\
\alpha_s>\theta_s-\theta_h >-\alpha_s\hspace{40px}The\; rate\;is\;[0,1].\hspace{186px}\\
\end{aligned}
$$
[Sun](https://www.geogebra.org/geometry/gpqfjudf)

## Single Scattering

单散射辐射是太阳射出的光线，(由于空气中分子或者气溶胶粒子; 排除地面的散射，单独计算)经过大气中某一点只发生一次散射的情况。下面将介绍如何计算、存储在LUT中，读取。  

> [地面散射](https://ebruneton.github.io/precomputed_atmospheric_scattering/atmosphere/functions.glsl.html#irradiance)

### Computation

假设太阳光在 q 点到达 p 点前，由于空气粒子发生散射。  
<font color=#888888>这里以 空气粒子 导致的 Rayleigh 散射为例，如果是气溶胶导致的 Mie 散射，其实是一样的，</font>  
<font color=#888888>把下述提到的 Rayleigh 换成 <del>我们米米</del> Mie 即可。</font>  
![ComputationPQ](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/singleScattering_Computation.png?raw=true)  
那么到达 p 点的辐射量 受到以下部分影响：

* 来自顶部大气的太阳辐射
* 太阳 和 q 点之间的 transmittance (即在大气顶部到达 q 点的太阳光的比例)
* q 点的 Rayleigh scattering coefficient (影响到达 q 点的光在任意方向上散射的比例)  
  可以理解为有多少散出去哩
* Rayleigh Phase Function (影响光在 q 点向 p 点散射多少比例)  
  可以理解为散出去的光有多少到达了 p 方向
* qp 之间的 transmittance (即在 q 点的散射光沿 qp 方向 到达 p 的比例)

因此，设有 指向太阳方向的单位向量 $w_s$ 即 $sun_{direction}$，以及：
$$
\begin{aligned}
&r=||\vec{op}||\\
&d=||\vec{pq}||\\
&\mu=\frac{\vec{op}\cdot\vec{pq}}{r*d}\\
&\mu_s=\frac{\vec{op}\cdot w_s}{r}\\
&\nu=\frac{\vec{pq}\cdot{w_s}}{d}\;这里\;\nu\;是\vec{pq}和太阳夹角,\,\\
&\nu是nu不是名扬天下的v\\
&r_d=||\vec{oq}||=\sqrt{d^2+2\mu *d*r+r^2}\\
&\mu_{s,d}=\frac{\vec{oq}\cdot w_s}{r_d}=
		\frac{(\vec{op}+\vec{pq})\cdot w_s}{r_d}=
		\frac{r*\mu_s+v*d}{r_d}
\end{aligned}
$$
Rayleigh Scattering 和 Mie Scattering 计算如下(这里先不带着 太阳辐照度和 Phase Function 计算，也暂时不考虑大气层底部，我们在后面进行处理)：
~~~C++
void ComputeSingleScatteringIntegrand(
	AtmosphereParameters atmosphere,
    TransmittanceTexture transmittance_texture,
    float r, float cosTheta, float cosTheta_s, float cosTheta_d, float d,
    bool ray_r_cosTheta_intersects_ground,
	out float3 rayleigh, out mie
)
    // cosTheta rayDir
    // cosTheta_s op zenth with sun
    // cosTheta_q pq zenth with sun
{
    float r_d = sqrt(d*d + 2*r*cosTheta*d + r*r);
    float cosTheta_s_d = (r*cosTheta_s + d*cosTheta_d) / r_d;
    // compute pq transmittance first
    // get sun transmittance on q
    // mul get sun_transmittance qp 
    float3 transmittance = 
        GetTransmittance(
    		atmosphere, transmittance_texture, r, cosTheta, d,
    		ray_r_cosTheta_intersects_ground) * 
        GetTransmittanceToSun(
    		atomsphere, transmittance_texture, r_d, cosTheta_s_d);
    rayleigh = transmittance * GetProfileDensity(
    	atmosphere.rayleigh_density, r_d - atmosphere.bottom_radius);
    mie = transmittance * GetProfileDensity(
    	atmosphere.mie_density, r_d - atmosphere.bottom_radius);
}
~~~

如果光线从给定方向 $\omega$，正好经过一次散射。这一次散射可能发生在 $p$ 和 射线$[p,\omega)$与大气边界交点 i 之间的任意一点 $q$。因此，$p$ 点 从 $\omega$ 方向的单散射 radiance总计 是 p 到 i 之间所有 <u>q 点到 p 点的单散射</u> 积分，所以我们得先知道 $||\vec{pi}||$
~~~C++
float DistanceToNearestAtmosphereBoundary(
	AtmosphereParameters atmosphere,
	float r, float cosTheta, bool ray_r_cosTheta_intersects_ground
)
{
    if(ray_r_cosTheta_intersects_ground)
    {
        return DistanceToBottomAtmosphereBoundary(atmosphere, r, cosTheta);
    }
    else
    {
        return DistanceToTopAtmosphereBoundary(atmosphere, r, cosTheta);
    }
}
~~~

单散射积分用 [梯形公式](https://en.wikipedia.org/wiki/Trapezoidal_rule)
~~~C++
void ComputeSingleScattering(
	AtmosphereParameters atmosphere,
    Transmittance transmittance_texture,
    float r, float cosTheta, float cosTheta_s, float cosTheta_d,
    bool ray_r_cosTheta_intersects_ground,
    out float3 rayleigh, out float3 mie
)
{
    assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
    assert(cosTheta >= -1.0 && cosTheta <= 1.0);
    assert(cosTheta_s >= -1.0 && cosTheta_s <= 1.0);
    assert(cosTheta_d >= -1.0 && cosTheta_d <= 1.0);

    const int SAMPLE_COUNT = 50;
    float dx = 
        DistanceToNearestAtmosphereBoundary(atmosphere, r, cosTheta,
			ray_r_cosTheta_intersects_ground) / (float)SAMPLE_COUNT;
    
    float3 rayleigh_sum = 0.0f;
    float3 mie_sum = 0.0f;
    for(int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        float d_i = (float)i * dx;
        
        float3 rayleigh_i;
        float3 mie_i;
        
        ComputeSingleScatteringIntegrand(atmosphere, transmittance_texture,
            r, cosTheta, cosTheta_s, cosTheta_d, d_i,
         	ray_r_cosTheta_intersects_ground, rayleigh_i, mie_i);
        
        float weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
    	rayleigh_sum += rayleigh_i * weight_i;
    	mie_sum += mie_i * weight_i;
    }
    rayleigh = rayleigh_sum * dx * atmosphere.solar_irradiance *
	atmosphere.rayleigh_scattering;
	mie = mie_sum * dx * atmosphere.solar_irradiance * 
		atmosphere.mie_scattering;
}
~~~

  
请注意，我们在`ComputeSingleScatteringIntegrand`中添加了太阳辐照度和散射系数项，这是我们之前忽略的部分，但没有添加 phase function - **phase function 项** 是在[渲染时](https://ebruneton.github.io/precomputed_atmospheric_scattering/atmosphere/functions.glsl.html#rendering)添加的，以提高角度精度。我们在这里提供它们是为了保证完整性：
~~~c++
float3 RayleighPhaseFunction(float cosTheta) {
  float3 k = 3.0 / (16.0 * PI * sr);
  return k * (1.0 + cosTheta * cosTheta);
}

float3 MiePhaseFunction(float g, float cosTheta) {
  float3 k = 3.0 / (8.0 * PI * sr) * (1.0 - g * g) / (2.0 + g * g);
  return k * (1.0 + cosTheta * cosTheta) / pow(1.0 + g * g - 2.0 * g * cosTheta, 1.5);
}
~~~

### Precomputation

`ComputeSingleScattering` 函数的计算开销很大，并且很多计算结果都是计算多重散射所需要的。因此我们可以把计算结果存储到一张 texture 里作为一张 LUT，所以我们需要一张 4D texture 来对应 4 个参数$(r,\mu,\mu_s,\nu)$。  
$$
\begin{aligned}
&r:distance\; of\;p\;from\;planet \; center.\\
&\mu:\cos\theta\;\;ray\; direction\; to\;zenith.\\
&\mu_s:\cos\theta_s\;\;sun\;direction\;to\;zenith.\\
&\nu:\cos\theta_{s_d}\;\; sun \; direction \; to \; ray\; direction.
\end{aligned}
$$
和 4D texture 对应 $(r,\mu,\mu_s,\nu)\rightarrow(u,v,w,z)$，下面的映射基于Eric 的 paper:[Precomputed Atmospheric Scattering](https://inria.hal.science/inria-00288758/en) 并有部分改进()：  
* $\mu:\cos\theta$ 的映射考虑到最近大气层的最小距离，将 $\mu$ 映射到 $[0,1]$，原始方法未能覆盖$[0,1]$所有值。  
* $\mu_s:\cos\theta_{sun}$ 的映射比 paper 里的要复杂，(原本映射是使用了为地球大气情况选择的临时参数)，该映射基于到达顶层大气边界的距离（适用于太阳光线），与 μ 映射类似，只使用一个临时（但可配置）的参数。并且，与原始定义一样，它在地平线附近提供了增加的采样率。
~~~C++
// (nu,mu_s,mu,r)
float4 GetScatteringTextureUvwzFromRMuMuSNu(
	AtmosphereParameters atmosphere,
	Length r, float mu, float mu_s, float nu,
	bool ray_r_mu_intersects_ground)
{
	assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	assert(mu >= -1.0 && mu <= 1.0);
	assert(mu_s >= -1.0 && mu_s <= 1.0);
	assert(nu >= -1.0 && nu <= 1.0);
	// 在地面的水平射线到大气层顶部的距离
	float H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
		atmosphere.bottom_radius * atmosphere.bottom_radius);

	// p 点到 planet 的切线长度
	float rho = max(sqrt(r*r - 
		atmosphere.bottom_radius * atmosphere.bottom_radius));
	float u_r = GetTextureCoordFromUnitRange(
		rho / H, SCATTERING_TEXTURE_R_SIZE );

	// 用平方公式判断射线(r,mu) 是否和地面相交，
	// 参照 RayIntersectsGround
	float r_mu = r * mu;
	float discriminant = 
		r_mu * r_mu - r * r + 
		atmosphere.bottom_radius * atmosphere.bottom_radius;
	float u_mu;
	if(ray_r_mu_intersects_ground)
	{
		// 射线到地表的距离，这里是指射线到达地表的距离
		// 其最小值和最大值覆盖(与地面相交)所有 mu 值，
		// 射线范围即：从 (r,-1) 到 (r,mu_horizon)
		float d = -r_mu - max(0.0f,sqrt(discriminant);
		float d_min = r - atmosphere.bottom_radius;
		float d_max = rho;
		// d_max == d_min 是 p 点在地面上
		u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(
			d_max == d_min ? 0.0 : (d - d_min)/(d_max - d_min),
			SCATTERING_TEXTURE_MU_SIZE / 2);
	}
	else
	{
		// 射线到大气层交点的距离，其最大最小值覆盖(不与地面相交)
		// 所有的 mu，
		// 射线范围即：从 (r,1) 到 (r,mu_horizon)
		float d = - r_mu + max(0.0f,sqrt(discriminant + H*H));
		float d_min = atmosphere.top_radius - r;
		float d_max = rho + H;
		u_mu = 0.5 + 0.5 * GetTextureCoordFromRange(
			(d-min)/(d_max-d_min), 
			SCATTERING_TEXTURE_MU_SIZE / 2);
	}
	float d = DistanceToTopAtmosphereBoundary(
		atmosphere, atmosphere.bottom_radius, mu_s);
	float d_min = atmosphere.top_radius - atmosphere.bottom_radius;
	float d_max = H;
	float a = (d-d_min)/(d_max-d_min);
	float D = DistanceToTopAtmosphereBoundary(
		atmosphere, atmosphere.bottom_radius, 
		atmosphere.mu_s_min);
	float A = (D-d_min)/(d_max-d_min);
	// 这是一种临时的函数，对于太阳天顶角 mu_s=mu_s_min 时等于 0
	// (因为此时 d=D 且 a=A)，
	// 对于 mu_s=1 时等于 1,(因为此时 d=d_min 且 a=0)，
	// 并在 mu_s​=0 附近具有较大的斜率，以获取地平线附近更多的纹理样本。
	float u_mu_s = GetTextureCoordFromUnitRange(
		max(1.0 - a/A, 0.0) / (1.0 + a),
		 SCATTERING_TEXTURE_MU_S_SIZE);
	float u_nu = (nu + 1.0) / 2.0;
	return vec4(u_nu, u_mu_s, u_mu, u_r);
}
~~~
[uvwz - GeoGebra](https://www.geogebra.org/geometry/hvfwsrc7)  
[Atmo_mu_s | Desmos](https://www.desmos.com/calculator/yyplzlmfvv?lang=zh-CN)  
在 $\mu_s=0$ 处斜率小，$\mu_s$ 会对应更多的 UV 值。  

逆变换如下：
~~~C++
void GetRMuMuSNuFromScatteringTextureUvwz(
	AtmosphereParameters atmosphere, float4 uvwz, 
	out float r, out float mu, out float mu_s, out float nu,
	out bool ray_r_mu_intersects_ground)
{
	// distance to topAtmosphere boundary for a horizontal ray at ground level.
	float H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
		atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to horizon
	// 从采样空间到数学 [0,1] 空间
	float rho = H * GetUnitRangeFromTextureCoord(uvwz.w,
		SCATTREING_TEXTURE_R_SIZE);
	r = sqrt(rho*rho + 
		atmosphere.bottom_radius * atmosphere.bottom_radius);

	if(uvwz.z < 0.5)
	{
		// Ray Cross Ground
		// 射线(r,mu)到地面的距离，并且其最大最小值覆盖所有射线和
		// 地面相交情况的 mu，从(r,-1)到(r,mu_horizon)
		// 可以得到 mu:
		float d_min = r - atmosphere.bottom_radius;
		float d_max = rho;
		float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
			1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE/2);
		mu = d == 0.0 ? -1.0 :
			clamp(-(rho*rho+d*d)/(2.0*r*d),-1.0,1.0);  
	    ray_r_mu_intersects_ground = true;
	}
	else
	{
		// 射线到大气层的距离，其最大最小值覆盖所有射线和大气层相交的
		// mu,从 (r,1) 到 (r,horizon)
		float d_min = atmosphere.top_radius - r;
		float d_max = rho + H;
		float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
			2.0*uvwz.z-1.0,SCATTERING_TEXTURE_MU_SIZE/2);
		mu = d == 0.0 ? 1.0 : 
			clamp( (H*H - rho*rho - d*d)/(2.0*r*d), -1.0,1.0 ); 
	    ray_r_mu_intersects_ground = false;
	}
	float x_mu_s = GetUnitRangeFromTextureCoord(uvwz.y,SCATTERING_TEXTURE_MU_S_SIZE);
	float d_min = atmsphere.top_radius - atmosphere.bottom_radius;
	float d_max = H;
	float D = DistanceToTopAtmosphereBoundary(
      atmosphere, atmosphere.bottom_radius, atmosphere.mu_s_min);
	float A = (D - d_min) / (d_max - d_min);
	float a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
	float d = d_min + min(a, A) * (d_max - d_min);
	mu_s = d == 0.0 ? 1.0 :
		clamp((H * H - d * d) / (2.0 * atmosphere.bottom_radius * d), -1.0, 1.0);

	// 注意这里 nu 不需要 Unit 和 texcoord 的变换，这和后面的 4D<->3D 相关
	nu = clamp(uvwz.x * 2.0 - 1.0, -1.0, 1.0);
}
~~~
这里为什么`d_min +  min(a, A)  *  (d_max - d_min);`，这个 `min(a, A)` 用以限制 $\mu_s$ 范围不超出 $\mu_{s_{min}}$
按照上述代码，我们需要一个 4D 纹理，但是实际上这不现实……因此我们要再做一次映射，在 3D 和 4D 之间。下述函数解释了如何从 3D 到 4D 纹理坐标$(r,\mu,\mu_s,\nu)$。这是通过 “unpacking” 把 x 分成两个纹理坐标来实现的。请注意在最后如何对 ν 参数进行了夹断。这是因为 ν 不是一个完全独立的变量：其取值范围取决于 μ 和 μs（通过从天顶、视角和太阳单位方向向量的笛卡尔坐标计算这一点是可以看出的），并且前面的函数暗示了这一点（如果不遵守这一约束，它们的断言可能会失败）。
~~~C++
void GetRMuMuSNuFromScatteringTextureFragCoord(
	AtmosphereParameters atmosphere, float3 frag_coord,
    out float r, out float mu, out float mu_s, out float nu,
    out bool ray_r_mu_intersects_ground)
{
	const float4 SCATTERING_TEXTURE_SIZE = float4(
		SCATTERING_TEXTURE_NU_SIZE - 1,
		SCATTERING_TEXTURE_MU_S_SIZE,
		SCATTERING_TEXTURE_MU_SIZE,
		SCATTERING_TEXTURE_R_SIZE);
	// 小于等于 的 最大整数值
	// 把 SCATTERING_TEXTURE_MU_S_SIZE 的 x 分成一组，共 若干组从 [0,SCATTERING_TEXTURE_NU_SIZE-1] 组，每组 SCATTERING_TEXTURE_MU_S_SIZE 个x对应同一个 nu
	float frag_coord_nu =
	    floor(frag_coord.x / float(SCATTERING_TEXTURE_MU_S_SIZE));
    // 对 x 取Mu_S余数
    // 把 x 变成多个 [0, SCATTERING_TEXTURE_MU_S_SIZE-1] 表示 mu_s 
	float frag_coord_mu_s =
	    mod(frag_coord.x, float(SCATTERING_TEXTURE_MU_S_SIZE));
	float4 uvwz =
	    float4(frag_coord_nu, frag_coord_mu_s, frag_coord.y, frag_coord.z) /
	        SCATTERING_TEXTURE_SIZE;
	GetRMuMuSNuFromScatteringTextureUvwz(
	    atmosphere, uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	// Clamp nu to its valid range of values, given mu and mu_s.
	nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
	    mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}
~~~
$$
\begin{aligned}
\cos(\alpha+\theta)=\cos\alpha\cdot\cos\theta-\sin\alpha\cdot\sin\theta\\
\cos(\alpha-\theta)=\cos\alpha\cdot\cos\theta+\sin\alpha\cdot\sin\theta
\end{aligned}
$$
<font color=#EF4F00>为什么 `SCATTRING_TEXTURE_NU_SIZE - 1`?   
这样 nu 的 coord 是 [0,1],  
但是 mu_s, mu, r 都是 [0,1)即$[0,1-\frac{1.0}{tex\_size}]$</font>  
<font color=#055000>如果要和后面`GetUnitRangeFromTextureCoord`对应，那在 0 处就错了，是 $\frac{-0.5}{tex\_size}$  
所以是否 frag_coord 的范围就是 $[\frac{0.5}{frag\_size},1-\frac{0.5}{frag\_size}]$，  
这样对 $r,\mu,\mu_s$ 的纹理变换都是正常范围和前面规定的 unit range <-> texcoord 一致，  
而 $\nu$ 因为使用 floor 函数，也是正常得到[0,nu_size-1]后面可以变为[0,1]。  
这样一切都对了</font>

<font color=red> **为什么** `SCATTRING_TEXTURE_NU_SIZE - 1` **?**  </font>
首先这里得到的 frag_coord_应该是[0,nu_size]，而 texcoord 即为 [0,1]
因为在 `GetRMuMuSNuFromScatteringTextureUvwz` 中，对 $\nu$ 直接进行 `uvwz.x * 2.0 - 1.0`
没做额外的映射所以最简单粗暴的[0,1]，
而且后面对 3D tex 采样有关 $\nu$ 的插值由我们手动进行，
不必担心因为 纹理采样的外推导致采样结果不准确，也不需要做 texcoord 和 unit range 两者之间的转换。

因为这里原文也没说 `frag_coord` 这个参数怎么处理，就……呃呃呃，我也懒得找了……  
但是不这么处理肯定错了啊

完整的单散射预计算方法如下：
~~~C++
void ComputeSingleScatteringTexture(AtmosphereParameters atmosphere,
	TransmittanceTexture transmittance_texture, float3 frag_coord,
	out float3 rayleigh, out float3 mie)
{
	float r;
	float mu;
	float mu_s;
	float nu;
	bool ray_r_mu_intersects_ground;
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
      r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	ComputeSingleScattering(atmosphere, transmittance_texture,
      r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
}
~~~
### LookUp
通过上述预计算的纹理的帮助，我们现在可以通过两次纹理查找获取点与最近的大气边界之间的散射效果（我们需要两个3D纹理查找来模拟单个4D纹理查找，并进行四线性插值；3D纹理坐标是使用在 `GetRMuMuSNuFromScatteringTextureFragCoord` 中定义的 3D-4D 映射的逆过程计算的）。
( 还有一点，因为我们前面所做的纹理映射可以保证在不同 $\nu$ 的纹理之间不会存在因为插值导致的错误）：
~~~C++
float3 GetScattering(AtmosphereParameters atmosphere, 
	AbstractScatteringTexture scattering_texture, float r, float mu, float mu_s, float nu,
	bool ray_r_mu_intersects_ground)
{
	float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(atmosphere,
		r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	// tex_coord_x 获取 真正采样所期望的 coord
	float tex_coord_x = uvwz.x * float(SCATTERING_TEXCOORD_NU_SIZE - 1);
	// 实际采样的起始点
	float tex_x = floor(tex_coord_x);
	// 两个采样做插值的lerp值
	float lerp = tex_coord_x - tex_x;
	// mu_s : nu 组内 mu_s 的结果所以要加 uvwz.y
	float3 uvw0 = float3( (tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE), 
		uvwz.z, uvwz.w);
	float3 uvw1 = float3( (tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE), 
		uvwz.z, uvwz.w);
	return float3(texture(scattering_texture, uvw0) * (1.0 - lerp) +
		texture(scattering_texture, uvw1) * lerp);
}
~~~
最后，在这里提供一个便捷的查找函数，它将在下一部分中很有用。该函数返回单次散射，其中包括相位函数，或者是 n 阶散射，其中 n > 1。它假设如果 `scattering_order` 严格大于1，那么 `multiple_scattering_texture` 对应于该散射阶数，其中包括了rayleigh 散射和 mie 散射，以及所有相位函数项。
~~~C++
float3 GetScattering(
    IN(AtmosphereParameters) atmosphere,
    IN(ReducedScatteringTexture) single_rayleigh_scattering_texture,
    IN(ReducedScatteringTexture) single_mie_scattering_texture,
    IN(ScatteringTexture) multiple_scattering_texture,
    float r, float mu, float mu_s, float nu,
    bool ray_r_mu_intersects_ground,
    int scattering_order) {
  if (scattering_order == 1) {
    float3 rayleigh = GetScattering(
        atmosphere, single_rayleigh_scattering_texture, r, mu, mu_s, nu,
        ray_r_mu_intersects_ground);
    float3 mie = GetScattering(
        atmosphere, single_mie_scattering_texture, r, mu, mu_s, nu,
        ray_r_mu_intersects_ground);
    return rayleigh * RayleighPhaseFunction(nu) +
        mie * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
  } else {
    return GetScattering(
        atmosphere, multiple_scattering_texture, r, mu, mu_s, nu,
        ray_r_mu_intersects_ground);
  }
}
~~~
## Multiple Scattering
多次散射辐射是指太阳光在大气中经过两次或更多次散射或从地面反射后到达某一点的光线。接下来的部分将描述我们如何计算它，如何将其存储在预先计算的纹理中以及如何读取它。

需要注意的是，与单次散射一样，我们在这里排除了最后一次反弹是在地面反射的光线路径。这些路径的贡献将在渲染时单独计算，以考虑实际地面反照率（对于预先计算的地面中间反射，我们使用平均均匀的反照率）。
### Computation
多次散射可以分解为双次散射、三次散射等的总和，其中每级对应于太阳光在大气中经过确切的2、3等次反弹后到达某一点的光。此外，每一级都可以从前一级计算出来。实际上，从方向 ω 经过n次反弹到达某一点 p 的光是对所有可能的最后一次反弹点 q 的积分，其中涉及从任意方向经过 n−1 次反弹到达 q 点的光。

这个描述表明，每个散射阶数都需要从前一个阶数计算一个三重积分<font color=#AFAFAF>（在方向 ω 中，沿着从 p 到最近的大气边界的线段上的所有点 q 进行 <u> 一次积分 </u>，以及在每个点 q 处对所有方向进行嵌套的<u>双重积分</u>）</font>。因此，如果我们想要从头计算每个阶数，我们将需要为双次散射进行三重积分，为三次散射进行六重积分(二次的三次+自己三次)，依此类推。这显然是低效的，因为存在大量冗余计算<font color=#AFAFAF> (阶数n的计算基本上会重新进行所有先前阶数的计算，导致总阶数的二次复杂度)</font>。相反，以下方法更为高效：
* 预计算单次散射存储在 texture 中
* 如果散射层级 n ≥ 2:
	* 用三重积分预先计算第 n 级的散射，其被积函数存在第 (n-1) 次散射纹理中。

这种策略避免了许多冗余计算，但并没有消除所有的冗余。例如，考虑下图中的点p和p′，以及计算在经过n次反弹后从方向ω到达这两点的光的必要计算。这些计算涉及到对在点q处沿着方向−ω散射的辐射L的评估，以及经过n−1次反弹后来自所有方向的辐射：  
![enter image description here](https://github.com/HollowEmiya/EmiyaPicGoRepo/blob/main/AtmosphereScattering/multipleScatteringCompute.png?raw=true)
因此，如果用如上的方法通过三重积分计算了第 n 级的散射，我们重复计算了 L 部分(对于q和最近大气边界方向−ω之间的所有点p)。
为了避免这种情况，从而提高多次散射计算的效率，我们将上述算法改进为: 
## 参考

[ebruneton/precomputed_atmospheric_scattering: This project provides a new implementation of our EGSR 2008 paper "Precomputed Atmospheric Scattering". (github.com)](https://github.com/ebruneton/precomputed_atmospheric_scattering)  
[ebruneton precomputed_atmospheric_scattering/](https://ebruneton.github.io/precomputed_atmospheric_scattering/)  
[Ebruneton PAS: Functions.glsl](https://ebruneton.github.io/precomputed_atmospheric_scattering/atmosphere/functions.glsl.html)

[预计算大气散射模型：原理与实现 - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/383020796) 

[^0]:我的地平线为什么紫了/蓝了: [大气散射：7.预计算算法之痛点 - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/365495154)

画数学图太好用了，GeoGebra 神！：[几何 - GeoGebra](https://www.geogebra.org/geometry)

[PicGo is Here | PicGo](https://picgo.github.io/PicGo-Doc/zh/guide/#picgo-is-here)
<!--stackedit_data:
eyJkaXNjdXNzaW9ucyI6eyJGUDZ1dU9HcGQ4Wno1NFdtIjp7In
N0YXJ0IjoyMzY2OSwiZW5kIjoyMzY5NSwidGV4dCI6InJheV9y
X211X2ludGVyc2VjdHNfZ3JvdW5kIn19LCJjb21tZW50cyI6ey
JKZjVSZ0JJeW5qVVBadTNIIjp7ImRpc2N1c3Npb25JZCI6IkZQ
NnV1T0dwZDhaejU0V20iLCJzdWIiOiJnaDo3MzQxOTk1NCIsIn
RleHQiOiLlsITnur/mmK/lkKblkozlnLDpnaLnm7jkuqQiLCJj
cmVhdGVkIjoxNzA2MTc4NjM0ODEzfX0sImhpc3RvcnkiOlstNz
IwNDU3OCw0ODgyODc2OTIsLTE2OTc1MzMwNjMsMTk1NjUxODg1
NCwtMTI4ODQ1MDQzNiwtNDk4MTc1NTgxLC0xMzMyMzI2NDAwLC
0zODY2NTk5NzIsLTIxMDA3MzI3MjksMjEyODQxNzQ3OSwtNDAw
OTI0NjkyLC0xNDY2ODk3OTMyLC0xNTY2NDcyMTksLTEyNzg4Nj
M2NTQsNTgwNTMyMjIxLC05NjkwNTYzNDMsMTEwOTM5ODk2MSwt
NTU4MTk5MTM0LC0xNDQ1NzAwNjI1LC0xNzY0MDU3NzMzXX0=
-->
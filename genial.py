import numpy as np

from enum import IntEnum
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import truncnorm

class Galaxy(object):    
    class Color(IntEnum):
        RED = 0
        GREEN = 1
        BLUE = 2
    class Kind(IntEnum):
        CENTRAL = 0
        SATELLITE = 1

class IntrinsicAlignment(object):
    
    def __init__(self, shape_params):
        # FIXME: Do this outside
        points = shape_params.index.levels
        shape = tuple([len(e) for e in points] + [len(shape_params.columns)])
        values = shape_params.values.reshape(shape)
        
        self._interp = RegularGridInterpolator(points, values)

    def axis_ratio(self, gal_color, redshift, mag):
        # values for clipping the pdf
        q_min, q_max = (0.001,1.0)
        r_min, r_max = (0.001,1.0)
        
        # interpolate shape parameters
        # FIXME: Do outside using xarray interpolation
        q_mean, r_mean, sigma = self._interp((gal_color, redshift, mag)).T
        
        # limits for clipping gaussian pdf    
        q_clip_lo = (q_min - q_mean) / sigma
        q_clip_hi = (q_max - q_mean) / sigma
        r_clip_lo = (r_min - r_mean) / sigma
        r_clip_hi = (r_max - r_mean) / sigma
        
        # draw parameters from pdf
        r = truncnorm.rvs(r_clip_lo, r_clip_hi, loc=r_mean, scale=sigma)
        q = truncnorm.rvs(q_clip_lo, q_clip_hi, loc=q_mean, scale=sigma)
        s = q*r

        return q, s

    def _random_vector(self, size):
        """Returns a randomly oriented unit vector"""
        phi = np.random.uniform(0, np.pi*2, size)

        cos_theta = np.random.uniform(-1, 1, size)
        theta = np.arccos(cos_theta)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        return np.array([x,y,z]).T

    def _random_perpendicular(self, v):
        """get random vector which is perpendicular to input vector"""
        vec_rnd = self._random_vector(len(v))
        vec_perp = np.cross(v, vec_rnd)

        vec_perp /= np.linalg.norm(vec_perp, axis=1)[None,:].T

        return vec_perp

    def orientations(self, gal_kind, gal_color, vg, vh, Ah, Ch, Jh):
        """recipe for pointing galaxies"""
        Ag = np.full(Ah.shape, np.nan)
        Cg = np.full(Ah.shape, np.nan)

        central_blue = (gal_kind == Galaxy.Kind.CENTRAL) & (gal_color == Galaxy.Color.BLUE)
        central_red  = (gal_kind == Galaxy.Kind.CENTRAL) & (gal_color == Galaxy.Color.RED)
        satellites   = (gal_kind == Galaxy.Kind.SATELLITE)

        # CENTRAL-BLUE (minor axis aligned with angular momentum vector, major axis random)
        Ag[central_blue] = self._random_perpendicular(Jh[central_blue])
        Cg[central_blue] = Jh[central_blue]

        # CENTRAL-RED (same major and minor axis as host halo)
        Ag[central_red] = Ah[central_red]
        Cg[central_red] = Ch[central_red]

        # SATELLITES (major axis pointing to halo center, minor axis lies on tangential plane)
        Ag[satellites] = vg[satellites] - vh[satellites] 
        
        satellites_Ag_zero = satellites & (np.linalg.norm(Ag, axis=1)[:,None]==0).T[0]

        Ag[satellites_Ag_zero==True] = self._random_vector(len(Ag[satellites_Ag_zero==True]))
        Ag[satellites_Ag_zero==False] /= np.linalg.norm(Ag[satellites_Ag_zero==False], axis=1)[:,None]


        Cg[satellites] = self._random_perpendicular(Ag[satellites])

        return Ag, Cg
    
    def _rotate(self, v, phi, theta):
        v = v.T

        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        return np.array([
             v[0]*cos_theta*cos_phi - v[1]*sin_phi + v[2]*sin_theta*cos_phi,
             v[0]*cos_theta*sin_phi + v[1]*cos_phi + v[2]*sin_theta*sin_phi,
            -v[0]*sin_theta                        + v[2]*cos_theta
        ]).T

    def _theta_Fisher(self, sigma):
        """
        random angle, following Fisher distribution
        ref: Joachimi et al. 2011, https://arxiv.org/abs/1305.5791

        """
        # numerically stable for low sigma, but unstable for sigma > 100000.
        # Stability might depend on mashine,
        # so we set random angles for sigma>1000, which is a very good approximation
        def theta_lo(sigma,r):
            s = sigma**(-2)
            tmp = 1. - (((r - 1.) / r) * np.exp(-2. * s))
            return np.arccos(1. + (np.log(r) + np.log(tmp)) / s)

        def theta_hi(r):
            return np.arccos(2*r-1)

        theta = np.zeros(sigma.size)
        r = np.random.uniform(0, 1, sigma.size)

        low = sigma < 1000
        hi = low==False

        theta[low] = theta_lo(sigma[low],r[low])
        theta[hi] = theta_hi(r[hi])
    
        return theta

    def _rand_3Dvec_Fisher(self, v, sigma):
        """# make randomized version of input vector"""
        
        # input vector
        vec = v / np.linalg.norm(v, axis=1)[None,:].T

        # angular position of input vector
        theta = np.arccos(vec[:, 2])
        phi = np.arctan2(vec[:, 1], vec[:, 0])

        # random angles
        theta_rand = self._theta_Fisher(sigma)
        phi_rand = np.random.uniform(0, 2*np.pi, size=sigma.size)

        # randomize vector along z-axis
        vz = np.repeat([[0,0,1]], sigma.size, axis=0)
        vz_rand = self._rotate(vz, phi_rand, theta_rand)

        # rotate vz_rand into direction of input vector
        vec_rand = self._rotate(vz_rand, phi, theta)

        return vec_rand

    def misalignment_parameters_v0(self, gal_kind, gal_color, vm, p_cr, p_cb, p_sr, p_sb):
        """ set width of Fisher distribution as function of galaxy properties """

        centrals_red  = (gal_kind == Galaxy.Kind.CENTRAL) & (gal_color == Galaxy.Color.RED)
        satellites_red   = (gal_kind == Galaxy.Kind.SATELLITE) & (gal_color == Galaxy.Color.RED)

        centrals_blue = (gal_kind == Galaxy.Kind.CENTRAL) & (gal_color == Galaxy.Color.BLUE)
        satellites_blue   = (gal_kind == Galaxy.Kind.SATELLITE) & (gal_color == Galaxy.Color.BLUE)

        sigma = np.full(vm.shape, np.nan)

        m0=-22.

        sigma[centrals_blue] = p_cb[0] + p_cb[1] * (vm[centrals_blue]/m0-1);
        sigma[satellites_blue] = p_sb[0] + p_sb[1] * (vm[satellites_blue]/m0-1);

        sigma[centrals_red] = p_cr[0] + p_cr[1] * (vm[centrals_red]/m0-1);
        sigma[satellites_red] = p_sr[0] + p_sr[1] * (vm[satellites_red]/m0-1);
                
        sigma[sigma<0.1]=0.1

        return sigma
    

    def randomize(self, Ag, Cg, sigma):
        """make randomized versions of input vectors A and C"""

        # randomize major amd minor axis
        Ar = self._rand_3Dvec_Fisher(Ag, sigma)
        Cr = self._rand_3Dvec_Fisher(Cg, sigma)

        # semi-minor axis: B_rand = (A_rand x C_rand)
        Br = np.cross(Ar, Cr)

        # normalize B_rand
        Br /= np.linalg.norm(Br, axis=1)[None, :].T

        # orthogonal minor axis C_rand = (B_rand x A_rand)
        Cr = np.cross(Br, Ar)

        return Ar, Cr
      
    def _e1e2e3(self, v):
        """tangential plane vectors (used for defining intrinsic shear)"""
        
        # normalized position vector
        pos_norm = v / np.linalg.norm(v, axis=1)[None, :].T

        pos_xy = pos_norm * np.array([1, 1, 0])
        pos_xy = pos_xy / np.linalg.norm(pos_xy, axis=1)[None, :].T
        pos_z = np.repeat([[0,0,1]], len(pos_norm), axis=0)

        # tangential plane vectors
        e1 = np.cross(pos_xy, pos_z)
        e2 = np.cross(e1, pos_norm)
        
        # line of sight vecots
        e3 = -pos_norm

        return e1, e2, e3
    
    def _project(self, v, e1, e2, e3) :
        """express input vector in base of tangential plane-los coord. system"""
      
        v_proj = np.array([
            v[:, 0]*e1[:, 0] + v[:, 1]*e1[:, 1] + v[:, 2]*e1[:, 2],
            v[:, 0]*e2[:, 0] + v[:, 1]*e2[:, 1] + v[:, 2]*e2[:, 2],
            v[:, 0]*e3[:, 0] + v[:, 1]*e3[:, 1] + v[:, 2]*e3[:, 2],
        ])
   
        return v_proj.T
  
    def ellipticity(self, v, q, s, A, C):
        """
        input:
          v: position vector
          q: |B|/|A|
          s: |C|/|A|
          A: major axis vector of 3D ellipsoid
          C: minor axis vector of 3D ellipsoid
        """

        # get tangential plane vectors
        #(ellipticity 'eps' is later defined with repsect to e1,e2)
        e1, e2, e3 = self._e1e2e3 (v)

        # project A, C on e1,e2,e3
        Ap = self._project(A, e1, e2, e3)
        Cp = self._project(C, e1, e2, e3)

        # normalize
        Ap /= np.linalg.norm(Ap, axis=1)[None,:].T
        Cp /= np.linalg.norm(Cp, axis=1)[None,:].T

        # semi-minor axis
        Bp = np.cross(Ap, Cp)

        # --- polarisation a la Joachimi 2014 ---
        #(arXiv:1305.5791 eq.2, see also arXiv:1203.6833, eq 3-5)
        Bp /= q[None,:].T
        Cp /= s[None,:].T

        S = np.swapaxes(np.array([Ap, Bp, Cp]), 0, 1)

        alpha_sq = S[:,0,2]**2 + S[:,1,2]**2 + S[:,2,2]**2

        k0 = S[:,0,2]*S[:,0,0] + S[:,1,2]*S[:,1,0] + S[:,2,2]*S[:,2,0]
        k1 = S[:,0,2]*S[:,0,1] + S[:,1,2]*S[:,1,1] + S[:,2,2]*S[:,2,1]

        W_inv_00 = S[:,0,0]*S[:,0,0] + S[:,1,0]*S[:,1,0] + S[:,2,0]*S[:,2,0] - k0*k0 / alpha_sq
        W_inv_01 = S[:,0,0]*S[:,0,1] + S[:,1,0]*S[:,1,1] + S[:,2,0]*S[:,2,1] - k0*k1 / alpha_sq
        W_inv_10 = S[:,0,1]*S[:,0,0] + S[:,1,1]*S[:,1,0] + S[:,2,1]*S[:,2,0] - k1*k0 / alpha_sq
        W_inv_11 = S[:,0,1]*S[:,0,1] + S[:,1,1]*S[:,1,1] + S[:,2,1]*S[:,2,1] - k1*k1 / alpha_sq

        # invert matrix
        det_W_inv = W_inv_00*W_inv_11 - W_inv_01*W_inv_10;

        W_00 =   W_inv_11 / det_W_inv;
        W_01 = - W_inv_10 / det_W_inv;
        W_10 = - W_inv_01 / det_W_inv;
        W_11 =   W_inv_00 / det_W_inv;
        
        # complex ellipticity
        det_W = 1./det_W_inv

        denom = W_00 + W_11 + 2*det_W**0.5

        eps1 = (W_00 - W_11) / denom
        eps2 = 2 * W_01 / denom

        return eps1, eps2

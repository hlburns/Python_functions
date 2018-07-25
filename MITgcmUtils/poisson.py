"""For computing rotational / divergent fluxes from MITgcm output"""

from pylab import *
import os.path
import warnings
from scipy import sparse

class solver:
    """Solves poisson problem using guass-seidel iteration"""

    def __init__(self, MITgcmmodel_parent, eps=40):
        
        self.eps = eps
        self.grid = MITgcmmodel_parent
        self.dx = self.grid.dxc
        self.dy = self.grid.dyc
        if (self.dx.ndim==1): # cartesian grid
            self.dx, self.dy = meshgrid(self.dx,self.dy)
        self.dx2 = (0.5*(roll(self.dx, -1, axis=1) + self.dx))**2
        self.dy2 = (0.5*(roll(self.dy, -1, axis=0) + self.dy))**2
        self.dx2dy2 = self.dx2 * self.dy2
        self.recip_dx2_dy2 = (self.dx2 + self.dy2)**-1            
        
    # pick a solver
    def solve(self, rhs):
        try:
            import pyamg
            self.use_pyamg=True
            self.solve_pyamg(rhs)
        except ImportError:
            warnings.warn('Coulding find pyamg, using slow python Jacobi solver')
            self.use_pyamg=False
            self.solve_jacobi(rhs)

    def solve_jacobi(self, rhs):
        """Solves the del^2 phi = rhs, subject to the boundary conditions
            grad phi = 0 at domain boundaries. Returns phi."""
        
        # needs to be padded on both sides, defined on p point
        self.phi = zeros((self.grid.Nx+2,self.grid.Ny+2))
        self.old_phi = self.phi.copy()
        self.rhs = rhs
        # doesn't seem to work properly
        self.normfac = sqrt(mean((0.5*self.recip_dx2_dy2*self.dx2dy2*self.rhs)**2))
        
        if isinstance(rhs, ma.masked_array):
            self.setup_masks(rhs.mask)
        else:
            self.setup_masks(zeros(rhs.shape), dtype('bool'))
        
        err = self.iterate()
        
        while err > self.eps:
            err = self.iterate()
            print err
    
    def compute_error(self):
        # how to normalize the error
        v = (self.phi - self.old_phi).flat
        return sqrt(dot(v,v)) / self.normfac
    
    def setup_masks(self,mask):
        # pad the mask
        self.mask = zeros((self.grid.Ny+2,self.grid.Nx+2), dtype('bool'))
        self.mask[1:-1,1:-1] = mask
        # make it wrap around
        self.mask = self.fix_periodicity(self.mask)
   
        # define four masks based on where the boundaries are located
        # western boundary
        if ~self.use_pyamg:
            self.eb = self.mask & ~roll(self.mask, -1, axis=1)
            self.wb = self.mask & ~roll(self.mask, 1, axis=1)
            self.sb = self.mask & ~roll(self.mask, -1, axis=0)
            self.nb = self.mask & ~roll(self.mask, 1, axis=0)

    def fix_periodicity(self, thevar=None):
        if thevar==None:
            thevar = self.phi
        thevar[0,:] = thevar[-2,:] # bottom
        thevar[-1,:] = thevar[1,:] # top
        thevar[:,0] = thevar[:,-2] # left
        thevar[:,-1] = thevar[:,1] # right
        return thevar

    def enforce_BCs(self):
        # make sure grad phi dot n = 0
        self.fix_periodicity()
        self.phi[self.wb] = self.phi[roll(self.wb, 1, axis=1)]
        self.phi[self.eb] = self.phi[roll(self.eb, -1, axis=1)]
        self.phi[self.nb] = self.phi[roll(self.nb, 1, axis=0)]
        self.phi[self.sb] = self.phi[roll(self.sb, -1, axis=0)]
            
    def iterate(self):
        # laplacian (centered diff):
        
        self.old_phi = self.phi.copy()
        
        self.phi[1:-1,1:-1] = 0.5 * self.recip_dx2_dy2 * (
                    + self.dy2*(self.phi[1:-1,0:-2] + self.phi[1:-1,2:])
                    + self.dx2*(self.phi[0:-2,1:-1] + self.phi[2:,1:-1])
                    - self.dx2dy2 * self.rhs )
        self.enforce_BCs()
                    
        return self.compute_error()
        
    def get_vectors(self):
        """Return vectors from phi, the potential function"""
        u = (self.phi[1:-1,1:-1] - self.phi[1:-1,0:-2]) / self.dx
        v = (self.phi[1:-1,1:-1] - self.phi[0:-2,1:-1]) / self.dy
        return ma.masked_array(u, self.mask[1:-1,1:-1]), ma.masked_array(v, self.mask[1:-1,1:-1])
        
    def solve_pyamg(self, rhs):
        import pyamg
        # construct the matrix form of the problem
        Nx = self.grid.Nx
        Ny = self.grid.Ny
        mask = rhs.mask
        
        self.phi = zeros((Ny+2,Nx+2))
        self.setup_masks(rhs.mask)
        
        if any(diff(self.grid.dxc)) | any(diff(self.grid.dyc)):
            warnings.warn('Detected unevenly spaced grid. solve_pyamg is not set up for this yet!')
        dx2 = self.grid.dxc[0]**2
        dy2 = self.grid.dyc[0]**2
        
        # these boundary indices are different from the ones above
        # they are wet points with land nearby
        nb = ~mask & roll(mask, -1, axis=0) # there is land directly north
        sb = ~mask & roll(mask, 1, axis=0)  # there is land directly south
        wb = ~mask & roll(mask, -1, axis=1) # there is land directly west
        eb = ~mask & roll(mask, 1, axis=1)  # there is land directly east
        N = Nx*Ny
        DX = pyamg.gallery.stencil_grid([[0,0,0],[1,-2,1],[0,0,0]], grid=(Ny,Nx), format='coo')
        DY = pyamg.gallery.stencil_grid([[0,1,0],[0,-2,0],[0,1,0]], grid=(Ny,Nx), format='coo')
        
        periodicEW = True
        periodicNS = False
        # left periodicity
        # (there are Ny overlapping X points at western boundary)
        row_j = Nx*arange(Ny)
        col_i = Nx*arange(Ny) + (Nx-1)
        olx_left = sparse.coo_matrix((ones(Ny), (row_j,col_i)), shape=(N,N))
        # right periodicity
        # (there are Ny overlapping X points at eastern boundary)
        row_j = Nx*arange(Ny) + (Nx-1)
        col_i = Nx*arange(Ny)
        olx_right = sparse.coo_matrix((ones(Ny), (row_j,col_i)), shape=(N,N))        
        # down periodicity
        # (there are Nx overlapping Y points at southern boundary)
        row_j = arange(Nx) + (N-Nx)
        col_i = arange(Nx) 
        oly_down = sparse.coo_matrix((ones(Nx), (row_j,col_i)), shape=(N,N))
        # right periodicity
        # (there are Nx overlapping Y points at northern boundary)
        row_j = arange(Nx) 
        col_i = arange(Nx) + (N-Nx)
        oly_up = sparse.coo_matrix((ones(Nx), (row_j,col_i)), shape=(N,N))

        if periodicEW:
            DXp = DX + olx_left + olx_right
        else:
            DXp = DX
        if periodicNS:
            DYp = DY + oly_up + oly_down
        else:
            DYp = DY
            
        # neumann BC at southern boundary
        s_bndry = sparse.coo_matrix((ones(Nx), (arange(Nx)+1,arange(Nx))), shape=(N,N))
        n_bndry = sparse.coo_matrix((ones(Nx), (arange(Nx)-1+(N-Nx),arange(Nx)+(N-Nx))), shape=(N,N))
        DYp = DYp + s_bndry + n_bndry



        # deal with boundaries in x
        dx_bndry_matrix = sparse.lil_matrix((N,N))
        for row in find(eb):
            DXp_row = DXp.getrow(row)
            col_right = Nx*floor(row/Nx)+mod(row+1,Nx)
            #col_right = row+1
            dx_bndry_matrix[row,DXp_row.indices] = -DXp_row.data
            dx_bndry_matrix[row,row] -= 1
            dx_bndry_matrix[row,col_right] += 1
        for row in find(wb):
            DXp_row = DXp.getrow(row)
            col_left = Nx*floor(row/Nx)+mod(row-1,Nx)
            dx_bndry_matrix[row,DXp_row.indices] = -DXp_row.data
            dx_bndry_matrix[row,row] -= 1
            dx_bndry_matrix[row,col_left] += 1
            
        # deal with boundaries in y
        dy_bndry_matrix = sparse.lil_matrix((N,N)) 
        for row in find(sb):
            DYp_row = DYp.getrow(row)
            col_up = floor(row/N) + mod(row+Nx,N)
            dy_bndry_matrix[row,DYp_row.indices] = -DYp_row.data
            dy_bndry_matrix[row,row] -= 1
            dy_bndry_matrix[row,col_up] += 1
        for row in find(nb):
            DYp_row = DYp.getrow(row)
            col_down = floor(row/N) + mod(row-Nx,N)
            dy_bndry_matrix[row,DYp_row.indices] = -DYp_row.data
            dy_bndry_matrix[row,row] -= 1
            dy_bndry_matrix[row,col_down] += 1
            
        #A = (DXp + dx_bndry_matrix)/dx2 + (DYp + dy_bndry_matrix)/dy2
        A = DXp/dx2 + DYp/dy2
        
        # don't operate on masked rows (fill with identity matrix)
        mask_matrix = sparse.lil_matrix((N,N)) # we will add rows to this as needed
        for row in find(mask):
            #print row
            Arow = A.getrow(row)
            mask_matrix[row,Arow.indices] = -Arow.data
            mask_matrix[row,row] = 1
        
        # ready to solve!
        self.A = A + mask_matrix
        self.ml = pyamg.smoothed_aggregation_solver(self.A, max_coarse=10)
        
        self.residuals=[]
        phi=self.ml.solve(b=rhs.flatten(), residuals=self.residuals,
            tol=1e-10, maxiter=500, accel='gmres')
        self.mynorm = norm( A.dot(phi) - rhs.flatten())

        phi.shape = (Ny,Nx)
        self.phi[1:-1,1:-1] = phi
        self.fix_periodicity() 
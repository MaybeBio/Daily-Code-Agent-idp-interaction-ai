#   This file is part of the markovmodel/deeptime repository.
#   Copyright (C) 2017, 2018 Computational Molecular Biology Group,
#   Freie Universitaet Berlin (GER)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import tensorflow as tf
import numpy as np
import scipy
import matplotlib.pyplot as plt

class VampnetTools(object):
    
    '''Wrapper for the functions used for the development of a VAMPnet.
    
    Parameters
    ----------
        
    epsilon: float, optional, default = 1e-10
        threshold for eigenvalues to be considered different from zero,
        used to prevent ill-conditioning problems during the inversion of the 
        auto-covariance matrices.
        
    k_eig: int, optional, default = 0
        the number of eigenvalues, or singular values, to be considered while
        calculating the VAMP score. If k_eig is higher than zero, only the top
        k_eig values will be considered, otherwise teh algorithms will use all
        the available singular/eigen values.
    '''

    def __init__(self, epsilon=1e-6, k_eig=0):
        self._epsilon = epsilon
        self._k_eig = k_eig
        
    @property
    def epsilon(self):
        return self._epsilon
   
    @epsilon.setter
    def epsilon(self, value):
        self._epsilon = value
    
    @property
    def k_eig(self):
        return self._k_eig
   
    @k_eig.setter
    def k_eig(self, value):
        self._k_eig = value


    def loss_VAMP(self, y_true, y_pred):
        '''Calculates the gradient of the VAMP-1 score calculated with respect
        to the network lobes. Using the shrinkage algorithm to guarantee that
        the auto-covariance matrices are really positive definite and that their
        inverse square-root exists. Can be used as a losseigval_inv_sqrt function
        for a keras model
        
        Parameters
        ----------
        y_true: tensorflow tensor.
            parameter not needed for the calculation, added to comply with Keras
            rules for loss fuctions format.
            
        y_pred: tensorflow tensor with shape [batch_size, 2 * output_size]
            output of the two lobes of the network
            
        Returns
        -------
        loss_score: tensorflow tensor with shape [batch_size, 2 * output_size].
            gradient of the VAMP-1 score
        '''
        
        # Remove the mean from the data
        x, y, batch_size, output_size = self._prep_data(y_pred)

        # Calculate the inverse root of the auto-covariance matrices, and the 
        # cross-covariance matrix
        matrices = self._build_vamp_matrices(x, y, batch_size)
        cov_00_ir, cov_11_ir, cov_01 = matrices
        
        vamp_matrix = tf.linalg.matmul(cov_00_ir, tf.linalg.matmul(cov_01, cov_11_ir))
        D, U, V = tf.linalg.svd(vamp_matrix, full_matrices=True, compute_uv=True)
        diag = tf.linalg.diag(D)
        
        # Base-changed covariance matrices
        x_base = tf.linalg.matmul(cov_00_ir, U)
        y_base = tf.linalg.matmul(cov_11_ir, V)

        # Calculate the gradients
        nabla_01 = tf.linalg.matmul(x_base, y_base, transpose_b=True)
        nabla_00 = -0.5 * tf.linalg.matmul(x_base, tf.linalg.matmul(diag, x_base,  transpose_b=True))
        nabla_11 = -0.5 * tf.linalg.matmul(y_base, tf.linalg.matmul(diag, y_base,  transpose_b=True))


        # Derivative for the output of both networks.
        x_der = 2 * tf.linalg.matmul(nabla_00, x) + tf.linalg.matmul(nabla_01, y)
        y_der = 2 * tf.linalg.matmul(nabla_11, y) + tf.linalg.matmul(nabla_01, x,  transpose_a=True)
        
        x_der = 1/(batch_size - 1) * x_der
        y_der = 1/(batch_size - 1) * y_der

        # Transpose back as the input y_pred was
        x_1d = tf.transpose(x_der)
        y_1d = tf.transpose(y_der)

        # Concatenate it again
        concat_derivatives = tf.concat([x_1d, y_1d], axis=-1)

        # Stops the gradient calculations of Tensorflow
        concat_derivatives = tf.stop_gradient(concat_derivatives)
        
        # With a minus because Tensorflow minimizes the loss-function
        loss_score = - concat_derivatives * y_pred

        return loss_score
    
    
    def loss_VAMP2_autograd(self, y_true, y_pred):
        '''Calculates the VAMP-2 score with respect to the network lobes. Same function
        as loss_VAMP2, but the gradient is computed automatically by tensorflow. Added
        after tensorflow 1.5 introduced gradients for eigenvalue decomposition and SVD
        
        Parameters
        ----------
        y_true: tensorflow tensor.
            parameter not needed for the calculation, added to comply with Keras
            rules for loss fuctions format.
            
        y_pred: tensorflow tensor with shape [batch_size, 2 * output_size]
            output of the two lobes of the network
            
        Returns
        -------
        loss_score: tensorflow tensor with shape [batch_size, 2 * output_size].
            gradient of the VAMP-2 score
        '''

        # Remove the mean from the data
        x, y, batch_size, output_size = self._prep_data(y_pred) 

        # Calculate the covariance matrices
        cov_01 = 1/(batch_size - 1) * tf.linalg.matmul(x, y, transpose_b=True)
        cov_00 = 1/(batch_size - 1) * tf.linalg.matmul(x, x, transpose_b=True) 
        cov_11 = 1/(batch_size - 1) * tf.linalg.matmul(y, y, transpose_b=True)

        # ✅ --- 关键修复：在这里添加正则化项以提高数值稳定性 ---
        # 获取矩阵的维度 n (output_size)
        n_dim = tf.shape(cov_00)[0]
        # 创建一个对角矩阵，其对角线元素是一个非常小的正数 (self.epsilon, 通常是 1e-7)
        regularization = tf.linalg.diag(tf.fill([n_dim], self.epsilon))
        
        # 将正则化项加到协方差矩阵的对角线上
        cov_00_reg = cov_00 + regularization
        cov_11_reg = cov_11 + regularization
        # --------------------------------------------------------

        # ✅ --- 现在，在稳定后的矩阵上执行求逆操作 ---
        # Calculate the inverse of the self-covariance matrices
        cov_00_inv = self._inv(cov_00_reg, ret_sqrt = True)
        cov_11_inv = self._inv(cov_11_reg, ret_sqrt = True) # <--- 这个调用现在是安全的了
        # ------------------------------------------------

        vamp_matrix = tf.linalg.matmul(tf.linalg.matmul(cov_00_inv, cov_01), cov_11_inv)

        vamp_score = tf.norm(vamp_matrix)

        return - tf.square(vamp_score)


    def loss_VAMP2(self, y_true, y_pred):
        '''Calculates the gradient of the VAMP-2 score calculated with respect
        to the network lobes. Using the shrinkage algorithm to guarantee that 
        the auto-covariance matrices are really positive definite and that their
        inverse square-root exists. Can be used as a loss function for a keras
        model
        
        Parameters
        ----------
        y_true: tensorflow tensor.
            parameter not needed for the calculation, added to comply with Keras
            rules for loss fuctions format.
            
        y_pred: tensorflow tensor with shape [batch_size, 2 * output_size]
            output of the two lobes of the network
            
        Returns
        -------
        loss_score: tensorflow tensor with shape [batch_size, 2 * output_size].
            gradient of the VAMP-2 score
        '''

        # Remove the mean from the data
        x, y, batch_size, output_size = self._prep_data(y_pred) 

        # Calculate the covariance matrices
        cov_01 = 1/(batch_size - 1) * tf.linalg.matmul(x, y, transpose_b=True)
        cov_10 = 1/(batch_size - 1) * tf.linalg.matmul(y, x, transpose_b=True)
        cov_00 = 1/(batch_size - 1) * tf.linalg.matmul(x, x, transpose_b=True) 
        cov_11 = 1/(batch_size - 1) * tf.linalg.matmul(y, y, transpose_b=True)
        
        # Calculate the inverse of the self-covariance matrices
        cov_00_inv = self._inv(cov_00)
        cov_11_inv = self._inv(cov_11)
        
        # Split the gradient computation in 2 parts for readability
        # These are reported as Eq. 10, 11 in the VAMPnets paper
        left_part_x = tf.linalg.matmul(cov_00_inv, tf.linalg.matmul(cov_01, cov_11_inv))
        left_part_y = tf.linalg.matmul(cov_11_inv, tf.linalg.matmul(cov_10, cov_00_inv))

        right_part_x = y - tf.linalg.matmul(cov_10, tf.linalg.matmul(cov_00_inv, x))
        right_part_y = x - tf.linalg.matmul(cov_01, tf.linalg.matmul(cov_11_inv, y))

        # Calculate the dot product of the two matrices
        x_der = 2/(batch_size - 1) * tf.linalg.matmul(left_part_x, right_part_x)
        y_der = 2/(batch_size - 1) * tf.linalg.matmul(left_part_y, right_part_y)

        # Transpose back as the input y_pred was
        x_1d = tf.transpose(x_der)
        y_1d = tf.transpose(y_der)

        # Concatenate it again
        concat_derivatives = tf.concat([x_1d,y_1d], axis=-1)

        # Stop the gradient calculations of Tensorflow
        concat_derivatives = tf.stop_gradient(concat_derivatives)

        # With a minus because Tensorflow maximizes the loss-function
        loss_score =  - concat_derivatives * y_pred
        
        return loss_score

    

    def metric_VAMP(self, y_true, y_pred):
        '''Returns the sum of the top k eigenvalues of the vamp matrix, with k
        determined by the wrapper parameter k_eig, and the vamp matrix defined
        as:
            V = cov_00 ^ -1/2 * cov_01 * cov_11 ^ -1/2
        Can be used as a metric function in model.fit()
        
        Parameters
        ----------
        y_true: tensorflow tensor.
            parameter not needed for the calculation, added to comply with Keras
            rules for loss fuctions format.
            
        y_pred: tensorflow tensor with shape [batch_size, 2 * output_size]
            output of the two lobes of the network
            
        Returns
        -------
        eig_sum: tensorflow float
            sum of the k highest eigenvalues in the vamp matrix
        '''

        # Remove the mean from the data
        x, y, batch_size, output_size = self._prep_data(y_pred)
        
        # Calculate the inverse root of the auto-covariance matrices, and the 
        # cross-covariance matrix
        matrices = self._build_vamp_matrices(x, y, batch_size)
        cov_00_ir, cov_11_ir, cov_01 = matrices
        
        # Calculate the VAMP matrix
        vamp_matrix = tf.linalg.matmul(cov_00_ir, tf.linalg.matmul(cov_01, cov_11_ir))
        
        # Select the K highest singular values of the VAMP matrix
        diag = tf.linalg.svd(vamp_matrix, compute_uv=False, full_matrices=False)
        cond = tf.greater(self.k_eig, 0)
        top_k_val = tf.nn.top_k(diag, k=self.k_eig)[0]
        
        # Sum the singular values
        eig_sum = tf.cond(cond, lambda: tf.reduce_sum(top_k_val), lambda: tf.reduce_sum(diag))
        
        return eig_sum


    def metric_VAMP2(self, y_true, y_pred):
        '''Returns the sum of the squared top k eigenvalues of the vamp matrix,
        with k determined by the wrapper parameter k_eig, and the vamp matrix
        defined as:
            V = cov_00 ^ -1/2 * cov_01 * cov_11 ^ -1/2
        Can be used as a metric function in model.fit()
        
        Parameters
        ----------
        y_true: tensorflow tensor.
            parameter not needed for the calculation, added to comply with Keras
            rules for loss fuctions format.
            
        y_pred: tensorflow tensor with shape [batch_size, 2 * output_size]
            output of the two lobes of the network
            
        Returns
        -------
        eig_sum_sq: tensorflow float
            sum of the squared k highest eigenvalues in the vamp matrix
        '''

        # Remove the mean from the data
        x, y, batch_size, output_size = self._prep_data(y_pred)
        
        # Calculate the inverse root of the auto-covariance matrices, and the 
        # cross-covariance matrix
        matrices = self._build_vamp_matrices(x, y, batch_size)
        cov_00_ir, cov_11_ir, cov_01 = matrices
        
        # Calculate the VAMP matrix
        vamp_matrix = tf.linalg.matmul(cov_00_ir, tf.linalg.matmul(cov_01, cov_11_ir))
        
        # Select the K highest singular values of the VAMP matrix
        diag = tf.linalg.svd(vamp_matrix, compute_uv=False, full_matrices=False)
        cond = tf.greater(self.k_eig, 0)
        top_k_val = tf.nn.top_k(diag, k=self.k_eig)[0]
        
        # Square the singular values and sum them
        pow2_topk = tf.reduce_sum(tf.multiply(top_k_val,top_k_val))
        pow2_diag = tf.reduce_sum(tf.multiply(diag,diag))
        eig_sum_sq = tf.cond(cond, lambda: pow2_topk, lambda: pow2_diag)
        
        return eig_sum_sq


    def estimate_koopman_op(self, trajs, tau):
        '''Estimates the koopman operator for a given trajectory at the lag time
            specified. The formula for the estimation is:
                K = C00 ^ -1 @ C01
                
        Parameters
        ----------
        traj: numpy array with size [traj_timesteps, traj_dimensions]
            or a list of trajectories
            Trajectory described by the returned koopman operator
            
        tau: int
            Time shift at which the koopman operator is estimated

        Returns
        -------
        koopman_op: numpy array with shape [traj_dimensions, traj_dimensions]
            Koopman operator estimated at timeshift tau

        '''
        if type(trajs) == list:
            traj = np.concatenate([t[:-tau] for t in trajs], axis = 0)
            traj_lag = np.concatenate([t[tau:] for t in trajs], axis = 0)
        else:
            traj = trajs[:-tau]
            traj_lag = trajs[tau:]
        
        c_0 = np.transpose(traj) @ traj
        c_tau = np.transpose(traj) @ traj_lag

        eigv, eigvec = np.linalg.eig(c_0)
        include = eigv > self._epsilon
        eigv = eigv[include]
        eigvec = eigvec[:,include]
        c0_inv = eigvec @ np.diag(1/eigv) @ np.transpose(eigvec)

        koopman_op = c0_inv @ c_tau
        return koopman_op


    def get_its(self, traj, lags):
        ''' Implied timescales from a trajectory estimated at a series of lag times.

        Parameters
        ----------
        traj: numpy array with size [traj_timesteps, traj_dimensions]
            trajectory data or a list of trajectories
            
        lags: numpy array with size [lag_times]
            series of lag times at which the implied timescales are estimated

        Returns
        -------
        its: numpy array with size [traj_dimensions - 1, lag_times]
            Implied timescales estimated for the trajectory.

        '''
        if type(traj) == list:
            outputsize = traj[0].shape[1]
        else:
            outputsize = traj.shape[1]
        its = np.zeros((outputsize-1, len(lags)))
 
        for t, tau_lag in enumerate(lags):
            koopman_op = self.estimate_koopman_op(traj, tau_lag)
            k_eigvals, k_eigvec = np.linalg.eig(np.real(koopman_op))
            k_eigvals = np.sort(np.absolute(k_eigvals))
            k_eigvals = k_eigvals[:-1]
            its[:,t] = (-tau_lag / np.log(k_eigvals))
            
        return its


    def get_ck_test(self, traj, steps, tau):
        ''' Chapman-Kolmogorov test for the koopman operator
        estimated for the given trajectory at the given lag times

        Parameters
        ----------
        traj: numpy array with size [traj_timesteps, traj_dimensions]
            trajectory data or a list of trajetories
            
        steps: int
            how many lag times the ck test will be evaluated at
            
        tau: int
            shift between consecutive lag times

        Returns
        -------
        predicted: numpy array with size [traj_dimensions, traj_dimensions, steps]
        estimated: numpy array with size [traj_dimensions, traj_dimensions, steps]
            The predicted and estimated transition probabilities at the
            indicated lag times

        '''

        if type(traj) == list:
            n_states = traj[0].shape[1]
        else:
            n_states = traj.shape[1]

        predicted = np.zeros((n_states, n_states, steps))
        estimated = np.zeros((n_states, n_states, steps))
 
        predicted[:,:,0] =  np.identity(n_states)
        estimated[:,:,0] =  np.identity(n_states)

        for vector, i  in zip(np.identity(n_states), range(n_states)):
            for n in range(1, steps):

                koop = self.estimate_koopman_op(traj, tau)

                koop_pred = np.linalg.matrix_power(koop,n)
                
                koop_est = self.estimate_koopman_op(traj, tau*n)
        
                predicted[i,:,n]= vector @ koop_pred
                estimated[i,:,n]= vector @ koop_est
        
              
        return [predicted, estimated]


    def estimate_koopman_constrained(self, traj, tau, th=0):
        ''' Calculate the transition matrix that minimizes the norm of the prediction
        error between the trajectory and the tau-shifted trajectory, using the
        estimate of the non-reversible koopman operator as a starting value.
        The constraints impose that all the values in the matrix are positive, and that
        the row sum equals 1. This is achieved using a COBYLA scipy minimizer. 

        Parameters
        ----------
        traj: numpy array with size [traj_timesteps, traj_dimensions]
            Trajectory described by the returned koopman operator
            
        tau: int
            Time shift at which the koopman operator is estimated
            
        th: float, optional, default = 0
            Parameter used to force the elements of the matrix to be higher than 0.
            Useful to prevent elements of the matrix to have small negative value
            due to numerical issues.

        Returns
        -------
        koop_positive: numpy array with shape [traj_dimensions, traj_dimensions]
            Koopman operator estimated at timeshift tau

        '''
        
        if type(traj) == list:
            raise Error('Multiple trajectories not supported')

        koop_init = self.estimate_koopman_op(traj, tau)
            
        n_states = traj.shape[1]
        
        
        rs = lambda k: np.reshape(k, (n_states, n_states))
        
        def errfun(k):
            diff_matrix = traj[tau:].T - rs(k) @ traj[:-tau].T
            return np.linalg.norm(diff_matrix)
        
        constr = []
        
        for n in range(n_states**2):
            # elements > 0
            constr.append({
                'type':'ineq',
                'fun': lambda x, n = n: x.flatten()[n] - th
                })
            # elements < 1
            constr.append({
                'type':'ineq',
                'fun': lambda x, n = n: 1 - x.flatten()[n] - th
                })
        
        for n in range(n_states):
            # row sum < 1
            constr.append({
                'type':'ineq',
                'fun': lambda x, n = n: 1 - np.sum(x.flatten()[n:n+n_states])
                })
            # row sum > 1
            constr.append({
                'type':'ineq',
                'fun': lambda x, n = n: np.sum(x.flatten()[n:n+n_states]) - 1
                })
        
        koop_positive = scipy.optimize.minimize(
            errfun, 
            koop_init,
            constraints = constr,
            method = 'COBYLA',
            tol = 1e-10,
            options = {'disp':False, 'maxiter':1e5},
            ).x
    
        return koop_positive


    def plot_its(self, its, lag, ylog=False):
        '''Plots the implied timescales calculated by the function
        'get_its'

        Parameters
        ----------
        its: numpy array
            the its array returned by the function get_its
        lag: numpy array
            lag times array used to estimate the implied timescales
        ylog: Boolean, optional, default = False
            if true, the plot will be a logarithmic plot, otherwise it
            will be a semilogy plot

        '''

        if ylog:
            plt.loglog(lag, its.T[:,::-1]);
            plt.loglog(lag, lag, 'k');
            plt.fill_between(lag, lag, 0.99, alpha=0.2, color='k');
        else:
            plt.semilogy(lag, its.T[:,::-1]);
            plt.semilogy(lag, lag, 'k');
            plt.fill_between(lag, lag, 0.99, alpha=0.2, color='k');
        plt.show()


    def plot_ck_test(self, pred, est, n_states, steps, tau):
        '''Plots the result of the Chapman-Kolmogorov test calculated by the function
        'get_ck_test'

        Parameters
        ----------
        pred: numpy array
        est: numpy array
            pred, est are the two arrays returned by the function get_ck_test
        n_states: int
        steps: int
        tau: int
            values used for the Chapman-Kolmogorov test as parameters in the function
            get_ck_test
        '''
        
        fig, ax = plt.subplots(n_states, n_states, sharex=True, sharey=True)
        for index_i in range(n_states):
            for index_j in range(n_states):
                
                ax[index_i][index_j].plot(range(0, steps*tau, tau),
                                          pred[index_i, index_j], color='b')
                
                ax[index_i][index_j].plot(range(0, steps*tau, tau),
                                          est[index_i, index_j], color='r', linestyle='--')
                
                ax[index_i][index_j].set_title(str(index_i+1)+ '->' +str(index_j+1),
                                               fontsize='small')
        
        ax[0][0].set_ylim((-0.1,1.1));
        ax[0][0].set_xlim((0, steps*tau));
        ax[0][0].axes.get_xaxis().set_ticks(np.round(np.linspace(0, steps*tau, 3)));
        plt.show()
        
    def _inv(self, x, ret_sqrt=False):
        '''Utility function that returns the inverse of a matrix, with the
        option to return the square root of the inverse matrix.
    
        Parameters
        ----------
        x: tensorflow tensor with shape [m,m]
            matrix to be inverted
            
        ret_sqrt: bool, optional, default = False
            if True, the square root of the inverse matrix is returned instead
    
        Returns
        -------
        x_inv: tensorflow tensor with shape [m,m]
            inverse of the original matrix
        '''
        # 1️⃣ 强制对称化（关键！）
        x = (x + tf.transpose(x)) / 2.0
    
        # 2️⃣ 添加正则化（关键！）
        n = tf.shape(x)[0]
        reg = self.epsilon * tf.eye(n, dtype=x.dtype)
        x = x + reg
    
        # 3️⃣ 特征分解
        eigval_all, eigvec_all = tf.linalg.eigh(x)
    
        # 4️⃣ 确保特征值非负（防数值误差）
        eigval_all = tf.maximum(eigval_all, self.epsilon)
    
        # 5️⃣ 构建逆或逆平方根
        if ret_sqrt:
            inv_sqrt_eigval = 1.0 / tf.sqrt(eigval_all)
            diag = tf.linalg.diag(inv_sqrt_eigval)
        else:
            inv_eigval = 1.0 / eigval_all
            diag = tf.linalg.diag(inv_eigval)
    
        # 6️⃣ 重构矩阵
        x_inv = tf.linalg.matmul(eigvec_all, tf.linalg.matmul(diag, eigvec_all, transpose_b=True))
        
        return x_inv


    def _prep_data(self, data):
        '''Utility function that transorms the input data from a tensorflow - 
        viable format to a structure used by the following functions in the
        pipeline.

        Parameters
        ----------
        data: tensorflow tensor with shape [b, 2*o]
            original format of the data

        Returns
        -------
        x: tensorflow tensor with shape [o, b]
            transposed, mean-free data corresponding to the left, lag-free lobe
            of the network
        
        y: tensorflow tensor with shape [o, b]
            transposed, mean-free data corresponding to the right, lagged lobe
            of the network
        
        b: tensorflow float32
            batch size of the data
        
        o: int
            output size of each lobe of the network
        
        '''

        shape = tf.shape(data)
        b = tf.cast(shape[0], tf.float32)
        o = shape[1]//2

        # Split the data of the two networks and transpose it
        x_biased = tf.transpose(data[:,:o])
        y_biased = tf.transpose(data[:,o:])

        # Subtract the mean
        x = x_biased - tf.reduce_mean(x_biased, axis=1, keepdims=True)
        y = y_biased - tf.reduce_mean(y_biased, axis=1, keepdims=True)

        return x, y, b, o
        

    def _build_vamp_matrices(self, x, y, b):
        '''Utility function that returns the matrices used to compute the VAMP
        scores and their gradients for non-reversible problems.

        Parameters
        ----------
        x: tensorflow tensor with shape [output_size, b]
            output of the left lobe of the network
            
        y: tensorflow tensor with shape [output_size, b]
            output of the right lobe of the network
            
        b: tensorflow float32
            batch size of the data
            
        Returns
        -------
        cov_00_inv_root: numpy array with shape [output_size, output_size]
            square root of the inverse of the auto-covariance matrix of x
            
        cov_11_inv_root: numpy array with shape [output_size, output_size]
            square root of the inverse of the auto-covariance matrix of y
        
        cov_01: numpy array with shape [output_size, output_size]
            cross-covariance matrix of x and y
            
        '''

        # Calculate the cross-covariance
        cov_01 = 1/(b-1) * tf.linalg.matmul(x, y, transpose_b=True)
        # Calculate the auto-correations
        cov_00 = 1/(b-1) * tf.linalg.matmul(x, x, transpose_b=True) 
        cov_11 = 1/(b-1) * tf.linalg.matmul(y, y, transpose_b=True)

        # Calculate the inverse root of the auto-covariance
        cov_00_inv_root = self._inv(cov_00, ret_sqrt=True)
        cov_11_inv_root = self._inv(cov_11, ret_sqrt=True)
 
        return cov_00_inv_root, cov_11_inv_root, cov_01

    
    def _build_vamp_matrices_rev(self, x, y, b):
        '''Utility function that returns the matrices used to compute the VAMP
        scores and their gradients for reversible problems. The matrices are
        transformed into symmetrical matrices by calculating the covariances
        using the mean of the auto- and cross-covariances, so that:
            cross_cov = 1/2*(cov_01 + cov_10)
        and:
            auto_cov = 1/2*(cov_00 + cov_11)


        Parameters
        ----------
        x: tensorflow tensor with shape [output_size, b]
            output of the left lobe of the network
            
        y: tensorflow tensor with shape [output_size, b]
            output of the right lobe of the network
            
        b: tensorflow float32
            batch size of the data

        Returns
        -------
        auto_cov_inv_root: numpy array with shape [output_size, output_size]
            square root of the inverse of the mean over the auto-covariance
            matrices of x and y
        
        cross_cov: numpy array with shape [output_size, output_size]
            mean of the cross-covariance matrices of x and y
        '''

        # Calculate the cross-covariances
        cov_01 = 1/(b-1) * tf.linalg.matmul(x, y, transpose_b=True)
        cov_10 = 1/(b-1) * tf.linalg.matmul(y, x, transpose_b=True)
        cross_cov = 1/2 * (cov_01 + cov_10) 
        # Calculate the auto-covariances
        cov_00 = 1/(b-1) * tf.linalg.matmul(x, x, transpose_b=True)
        cov_11 = 1/(b-1) * tf.linalg.matmul(y, y, transpose_b=True) 
        auto_cov = 1/2 * (cov_00 + cov_11) 

        # Calculate the inverse root of the auto-covariance
        auto_cov_inv_root = self._inv(auto_cov, ret_sqrt=True)
 
        return auto_cov_inv_root, cross_cov



    #### EXPERIMENTAL FUNCTIONS ####


    def _loss_VAMP_sym(self, y_true, y_pred):
        '''WORK IN PROGRESS        
        
        Calculates the gradient of the VAMP-1 score calculated with respect
        to the network lobes. Using the shrinkage algorithm to guarantee that 
        the auto-covariance matrices are really positive definite and that their
        inverse square-root exists. Can be used as a loss function for a keras
        model. The difference with the main loss_VAMP function is that here the
        matrices C00, C01, C11 are 'mixed' together:
        
        C00' = C11' = (C00+C11)/2
        C01 = C10 = (C01 + C10)/2
        
        There is no mathematical reasoning behind this experimental loss function.
        It performs worse than VAMP-2 with regard to the identification of processes,
        but it also helps the network to converge to a transformation that separates
        more neatly the different states
        
        Parameters
        ----------
        y_true: tensorflow tensor.
            parameter not needed for the calculation, added to comply with Keras
            rules for loss fuctions format.
            
        y_pred: tensorflow tensor with shape [batch_size, 2 * output_size]
            output of the two lobes of the network
            
        Returns
        -------
        loss_score: tensorflow tensor with shape [batch_size, 2 * output_size].
            gradient of the VAMP-1 score
        '''

        # Remove the mean from the data
        x, y, batch_size, output_size = self._prep_data(y_pred)
        
        # Calculate the inverse root of the auto-covariance matrix, and the 
        # cross-covariance matrix
        cov_00_ir, cov_01  = self._build_vamp_matrices_rev(x, y, batch_size)

        vamp_matrix = tf.linalg.matmul(cov_00_ir, tf.linalg.matmul(cov_01, cov_00_ir))

        D, U, V = tf.linalg.svd(vamp_matrix, full_matrices=True, compute_uv=True)
        diag = tf.linalg.diag(D)
        
        # Base-changed covariance matrices
        x_base = tf.linalg.matmul(cov_00_ir, U)
        y_base = tf.linalg.matmul(V, cov_00_ir, transpose_a=True)

        # Derivative for the output of both networks.
        nabla_01 = tf.linalg.matmul(x_base, y_base)
        nabla_00 = -0.5 * tf.linalg.matmul(x_base, tf.linalg.matmul(diag, x_base, transpose_b=True))

        # Derivative for the output of both networks.
        x_der = 2/(batch_size - 1) * (tf.linalg.matmul(nabla_00, x) + tf.linalg.matmul(nabla_01, y))
        y_der = 2/(batch_size - 1) * (tf.linalg.matmul(nabla_00, y) + tf.linalg.matmul(nabla_01, x))

        # Transpose back as the input y_pred was
        x_1d = tf.transpose(x_der)
        y_1d = tf.transpose(y_der)

        # Concatenate it again
        concat_derivatives = tf.concat([x_1d,y_1d], axis=-1)

        # Stop the gradient calculations of Tensorflow
        concat_derivatives = tf.stop_gradient(concat_derivatives)

        # With a minus because Tensorflow maximizes the loss-function
        loss_score = - concat_derivatives * y_pred
        
        return loss_score



    def _metric_VAMP_sym(self, y_true, y_pred):
        '''Metric function relative to the _loss_VAMP_sym function.
        
        Parameters
        ----------
        y_true: tensorflow tensor.
            parameter not needed for the calculation, added to comply with Keras
            rules for loss fuctions format.
            
        y_pred: tensorflow tensor with shape [batch_size, 2 * output_size]
            output of the two lobes of the network
            
        Returns
        -------
        eig_sum: tensorflow float
            sum of the k highest eigenvalues in the vamp matrix
        '''

        # Remove the mean from the data
        x, y, batch_size, output_size = self._prep_data(y_pred)
        
        # Calculate the inverse root of the auto-covariance matrices, and the 
        # cross-covariance matrix
        cov_00_ir, cov_01  = self._build_vamp_matrices_rev(x, y, batch_size)
        
        # Calculate the VAMP matrix
        vamp_matrix = tf.linalg.matmul(cov_00_ir, tf.linalg.matmul(cov_01, cov_00_ir))
        
        # Select the K highest singular values of the VAMP matrix
        diag = tf.linalg.svd(vamp_matrix, compute_uv=False, full_matrices=False)
        cond = tf.greater(self.k_eig, 0)
        top_k_val = tf.nn.top_k(diag, k=self.k_eig)[0]

        # Sum the singular values
        eig_sum = tf.cond(cond, lambda: tf.reduce_sum(top_k_val), lambda: tf.reduce_sum(diag))
        
        return eig_sum
    


    def _estimate_koopman_op(self, traj, tau):
        '''Estimates the koopman operator for a given trajectory at the lag time
            specified. The formula for the estimation is:
                K = C00 ^ -1/2 @ C01 @ C11 ^ -1/2
                
        Parameters
        ----------
        traj: numpy array with size [traj_timesteps, traj_dimensions]
            Trajectory described by the returned koopman operator
            
        tau: int
            Time shift at which the koopman operator is estimated

        Returns
        -------
        koopman_op: numpy array with shape [traj_dimensions, traj_dimensions]
            Koopman operator estimated at timeshift tau

        '''

        c_0 = traj[:-tau].T @ traj[:-tau]
        c_1 = traj[tau:].T @ traj[tau:]
        c_tau = traj[:-tau].T @ traj[tau:]

        eigv0, eigvec0 = np.linalg.eig(c_0)
        include0 = eigv0 > self._epsilon
        eigv0_root = np.sqrt(eigv0[include0])
        eigvec0 = eigvec0[:,include0]
        c0_inv_root = eigvec0 @ np.diag(1/eigv0_root) @ eigvec0.T

        eigv1, eigvec1 = np.linalg.eig(c_1)
        include1 = eigv1 > self._epsilon
        eigv1_root = np.sqrt(eigv1[include1])
        eigvec1 = eigvec1[:,include1]
        c1_inv_root = eigvec1 @ np.diag(1/eigv1_root) @ eigvec1.T

        koopman_op = c0_inv_root @ c_tau @ c1_inv_root
        return koopman_op

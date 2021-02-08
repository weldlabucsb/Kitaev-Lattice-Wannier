function [] = wannier_plots(wannier_func_realspace)
%WANNIER_PLOTS plot the real space wannier functions
%   Just to clean up the main Wannier Function
    fontsize = 26;
    figure
    surf(real(wannier_func_realspace));
    xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    ylabel('Y Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    zlab = ['Re(Wannier Func)'];
    zlabel(zlab, 'interpreter','latex','fontsize',fontsize);
    figure
    surf(imag(wannier_func_realspace));
    xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    ylabel('Y Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    zlab = ['Im(Wannier Func)'];
    zlabel(zlab, 'interpreter','latex','fontsize',fontsize);
    figure
    surf(abs(wannier_func_realspace));
    set(gca,'zscale','log');
%     surf(X,Y,abs(wannier_func_realspace));
    xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    ylabel('Y Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    zlab = ['Abs(Wannier)'];
    title('Wannier Func 1','interpreter','latex','fontsize',fontsize);
    zlabel(zlab, 'interpreter','latex','fontsize',fontsize);
    figure
    surf(angle(wannier_func_realspace));
    xlabel('X Pos., [$\lambda_l$]','interpreter','latex');
    ylabel('Y Pos., [$\lambda_l$]','interpreter','latex');
    zlab = ['Arg(Wannier Func)'];
    zlabel(zlab, 'interpreter','latex');
end


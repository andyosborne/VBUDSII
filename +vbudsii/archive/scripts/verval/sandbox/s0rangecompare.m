load fast_s0

subplot(1,2,1)
loglog(p.fineGroupDef(1:end-1), [rm1.Cell(1).spectralFlux r10.Cell(1).spectralFlux])

xlabel('energy (eV)')
ylabel('flux')
legend('s0 = -1','s0 = 10')

subplot(1,2,2)
loglog(p.fineGroupDef(1:end-1), [rm1.Cell(2).spectralFlux r10.Cell(2).spectralFlux])
xlabel('energy (eV)')
ylabel('flux')
legend('s0 = -1','s0 = 10')

function check_tt()
if (exist('tt_tensor', 'file')==0)
    if (exist('TT-Toolbox-small', 'dir')==0)
        if (exist('TT-Toolbox-small.zip', 'file')==0)
            try
                fprintf('TT-Toolbox is not found, downloading...\n');
                opts = weboptions; opts.CertificateFilename=(''); 
                websave('TT-Toolbox-small.zip', 'https://github.com/oseledets/TT-Toolbox/archive/small.zip', opts);
            catch ME
                warning(ME.identifier, '%s. Automatic download from github failed. Trying my website...', ME.message);
                try
                    opts = weboptions; opts.CertificateFilename=(''); 
                    websave('TT-Toolbox-small.zip', 'https://people.bath.ac.uk/sd901/als-cross-algorithm/ttsmall.zip', opts);
                catch ME2
                    error('%s. Automatic download failed. Please download TT-Toolbox from https://github.com/oseledets/TT-Toolbox', ME2.message);
                end
                fprintf('Success!\n');
            end
        end
        try
            unzip('TT-Toolbox-small.zip');
        catch ME
            error('%s. Automatic unzipping failed. Please extract TT-Toolbox-small.zip here', ME.message);
        end
    end
    cd('TT-Toolbox-small');
    setup;
    cd('..');
end

if (exist('tamen', 'file')==0)
    if (exist('tamen-master', 'dir')==0)
        if (exist('tamen-master.zip', 'file')==0)
            try
                fprintf('tamen is not found, downloading...\n');
                opts = weboptions; opts.CertificateFilename=(''); 
                websave('tamen-master.zip', 'https://github.com/dolgov/tamen/archive/refs/heads/master.zip', opts);
            catch ME
                error('%s. Automatic download failed. Please download the tAMEn package from https://github.com/dolgov/tamen', ME.message);
            end
        end
        try
            unzip('tamen-master.zip');
        catch ME
            error('%s. Automatic unzipping failed. Please extract tamen-master.zip here', ME.message);
        end
    end
    cd('tamen-master');
    addpath(pwd);
    cd('..');
end
end

function identity_replay(base)
% CUSTODIAL EXACT REPRESENTATION AUDIT. Never call the author's example
% script, MI routines, plotters, or classifiers. No response array, mean,
% standard deviation, difference, correlation, or ranking is printed/written.
% All choices are fixed in identity_check_plan.json before this function runs.

  addpath(fullfile(base, 'provenance', 'source_code'));
  addpath(fullfile(base, 'provenance', 'source_code_lf'));
  raw = jsondecode(fileread(fullfile(base, 'sealed', 'raw_mixed', 'fig1_sfp1_replicates.json')));
  targets = jsondecode(fileread(fullfile(base, 'sealed', 'raw_mixed', 'fig2_stress_type_replicates.json')));
  targets = targets.sfp1.gluc;
  labels = {'deposited_nucLoc_source_normalizer_only', ...
            'recomputed_ratio_source_normalizer_only', ...
            'deposited_nucLoc_worked_example', ...
            'recomputed_ratio_worked_example'};
  cases = cell(1, 48);
  k = 0;
  for r = 1:6
    rep = sprintf('rep%d', r);
    experiment = raw.(rep);
    origin = experiment.general.origin;
    times = experiment.general.times;
    representations = {experiment.GFP.nucLoc, ...
                       experiment.GFP.max5 ./ experiment.GFP.median};
    for v = 1:4
      original = representations{1 + mod(v - 1, 2)};
      succeeded = true;
      processed = [];
      try
        if v <= 2
          processed = normMeanOffset(original, origin, 20, 2);
        else
          % Exact source example_script.m lines 40-62, not a new fitted law.
          midtimes = median(times);
          synced = NaN(size(original));
          for c = 1:size(original, 1)
            mask = ~isnan(original(c, :));
            synced(c, :) = interp1(times(c, mask), original(c, mask), midtimes, 'linear');
          end
          [processed, ignored_cell_count] = filterCellsByNucLoc(synced, 0.5, 1, 0.8, origin);
          clear ignored_cell_count;
        end
      catch
        % Do not return data-dependent exception payloads or partial arrays.
        succeeded = false;
      end
      for t = 1:2
        target_rep = sprintf('rep%d', t);
        target = targets.(target_rep).GFP;
        k = k + 1;
        item = struct();
        item.source_group_id = ['1b38366b-7473-4ddf-9cf9-26e625533686:/' rep];
        item.target_group_id = ['721ef5f2-b858-4257-8cc8-5eea9a8599bb:/sfp1/gluc/' target_rep];
        item.variant_id = labels{v};
        item.known_duplicate_control = (r == 6 && t == 1);
        item.replay_completed = succeeded;
        item.same_central_window_shape = false;
        item.exact_double_equality = false;
        item.exact_source_decimal_token_equality = false;
        if succeeded && origin > 20 && target.origin > 20 && ...
            origin + 19 <= size(processed, 2) && ...
            target.origin + 19 <= size(target.nucLocNorm, 2)
          lhs = processed(:, origin - 20:origin + 19);
          rhs = target.nucLocNorm(:, target.origin - 20:target.origin + 19);
          item.same_central_window_shape = isequal(size(lhs), size(rhs));
          if item.same_central_window_shape && any(isfinite(rhs(:)))
            item.exact_double_equality = isequaln(lhs, rhs);
            % JSONlab-1.5 savejson.m: fixed default FloatFormat='%.10g'.
            % No adjustable numerical tolerance, row permutation or offset.
            item.exact_source_decimal_token_equality = strcmp( ...
                sprintf('%.10g,', lhs'), sprintf('%.10g,', rhs'));
          end
        end
        cases{k} = item;
      end
      clear original processed synced lhs rhs;
    end
  end
  report = struct();
  report.schema_version = 1;
  report.runtime = ['GNU Octave ' version()];
  report.response_values_or_statistics_released = false;
  report.native_models_fitted_or_scored = false;
  report.response_projections_written = false;
  report.cases = cases;
  % The only output operation in this custodian replay emits metadata only.
  fprintf('%s\n', jsonencode(report));
end

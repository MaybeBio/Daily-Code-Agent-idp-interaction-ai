from experiment_helpers import t as c
from train_jwarm import *
engine = c

def fit_stage(model, forward, splits, directory, args, device, *, stage, epochs, patience, lr, wd, stage_signature):
    directory.mkdir(parents=True, exist_ok=True)
    best_path, last_path = (directory / 'best.pt', directory / 'last.pt')
    if last_path.exists() and (not args.resume):
        raise RuntimeError(f'Existing stage {directory}; use --resume or a new output root')
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=3)
    start, stale, history, elapsed = (0, 0, [], 0.0)
    if last_path.exists():
        ckpt = load_pt(last_path)
        if ckpt['signature'] != stage_signature:
            raise RuntimeError(f'Stage configuration differs: {directory}')
        if ckpt['complete']:
            model.load_state_dict(ckpt['best_state'])
            if not best_path.exists():
                save_pt(best_path, dict(state=ckpt['best_state'], best_epoch=ckpt['best_epoch'], metrics=ckpt['best_metrics'], signature=stage_signature))
            logging.info('Reuse J_WARM %s best epoch=%d AP=%.6f', stage, ckpt['best_epoch'], ckpt['best_metrics']['ap'])
            return ckpt
        model.load_state_dict(ckpt['model'])
        optimizer.load_state_dict(ckpt['optimizer'])
        scheduler.load_state_dict(ckpt['scheduler'])
        restore_rng(ckpt['rng'], device)
        start, stale, history, elapsed = (ckpt['epoch'], ckpt['stale'], ckpt['history'], ckpt['seconds'])
        best_state, best_metric, best_epoch = (ckpt['best_state'], ckpt['best_metrics'], ckpt['best_epoch'])
    else:
        z = predict(model, forward, splits['val'], args.batch_size)
        best_metric = metrics(splits['val'].labels, z)
        best_state, best_epoch = (cpu_tree(model.state_dict()), 0)
        history.append(dict(epoch=0, train_bce=None, val_ap=best_metric['ap'], val_auc=best_metric['auc'], seconds=0.0))
    train = splits['train']
    for epoch in range(start + 1, epochs + 1):
        clock = time.monotonic()
        model.train()
        loss_total = 0.0
        for step, take in enumerate(batches(len(train), args.batch_size, args.seed + epoch)):
            optimizer.zero_grad(set_to_none=True)
            z = forward(model, 'train', take)
            y = torch.tensor(train.labels[take], dtype=torch.float32, device=device)
            loss = F.binary_cross_entropy_with_logits(z, y)
            if not torch.isfinite(loss):
                raise RuntimeError(f'{stage} loss is NaN/Inf')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
            optimizer.step()
            loss_total += float(loss.detach().cpu()) * len(take)
            if step and step % 200 == 0:
                logging.info('%s epoch=%d step=%d/%d BCE=%.5f', stage, epoch, step, math.ceil(len(train) / args.batch_size), loss_total / min((step + 1) * args.batch_size, len(train)))
        z = predict(model, forward, splits['val'], args.batch_size)
        val = metrics(splits['val'].labels, z)
        train_z = predict(model, forward, train, args.batch_size)
        train_eval = metrics(train.labels, train_z)
        train_eval_bce = float(np.mean(np.logaddexp(0, train_z) - train.labels * train_z))
        scheduler.step(val['ap'])
        if val['ap'] > best_metric['ap'] + 1e-06:
            best_state, best_metric, best_epoch, stale = (cpu_tree(model.state_dict()), val, epoch, 0)
        else:
            stale += 1
        seconds = time.monotonic() - clock
        elapsed += seconds
        history.append(dict(epoch=epoch, train_bce=loss_total / len(train), val_ap=val['ap'], val_auc=val['auc'], val_brier=val['brier'], train_eval_ap=train_eval['ap'], train_eval_bce=train_eval_bce, val_eval_bce=float(np.mean(np.logaddexp(0, z) - splits['val'].labels * z)), seconds=seconds))
        complete = stale >= patience or epoch == epochs
        ckpt = dict(signature=stage_signature, complete=complete, epoch=epoch, stale=stale, model=cpu_tree(model.state_dict()), optimizer=cpu_tree(optimizer.state_dict()), scheduler=scheduler.state_dict(), rng=rng_state(device), best_state=best_state, best_metrics=best_metric, best_epoch=best_epoch, history=history, seconds=elapsed, device=str(device))
        save_pt(last_path, ckpt)
        save_pt(best_path, dict(state=best_state, best_epoch=best_epoch, metrics=best_metric, signature=stage_signature))
        pd.DataFrame(history).to_csv(directory / 'history.tsv', sep='\t', index=False)
        logging.info('%s epoch=%d/%d BCE=%.5f valAP=%.6f valAUC=%.6f best=%d seconds=%.1f', stage, epoch, epochs, loss_total / len(train), val['ap'], val['auc'], best_epoch, seconds)
        release(device)
        if complete:
            break
    model.load_state_dict(best_state)
    return ckpt

def fit_joint(a, inputs, out, signature, device):
    t = inputs.training
    epochs = t['profile_epochs']
    patience = t['residual_patience']
    model = inputs.model(a, device, warm=True)
    require(torch.count_nonzero(model.profile.weights) == 0, 'W must initialize at zero')
    optimizer = torch.optim.AdamW(model.parameter_groups(t))
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=3)
    directory = out / 'joint'
    directory.mkdir(parents=True, exist_ok=True)
    last = directory / 'last.pt'
    best_path = directory / 'best.pt'
    c.seed_all(a.seed + 200)
    start = stale = 0
    elapsed = 0.0
    history = []
    if last.exists():
        require(a.resume, 'Use --resume to continue')
        ck = c.load_pt(last)
        require(ck['signature'] == signature, 'Training signature changed')
        if ck['complete']:
            model.load_state_dict(ck['best_state'])
            return (model, ck)
        model.load_state_dict(ck['model'])
        optimizer.load_state_dict(ck['optimizer'])
        scheduler.load_state_dict(ck['scheduler'])
        engine.restore_rng(ck['rng'], device)
        start, stale, elapsed, history = (ck['epoch'], ck['stale'], ck['seconds'], ck['history'])
        best_state, best, best_epoch = (ck['best_state'], ck['best_metrics'], ck['best_epoch'])
    else:
        take = np.arange(min(4, len(inputs.splits['val'])))
        model.eval()
        with torch.no_grad():
            z = model(*inputs.batch('val', take, device)).cpu().numpy()
        require(np.allclose(z, inputs.explicit['val'][take], atol=1e-05, rtol=1e-05), 'Warm start differs from pretrained explicit model')
        val = engine.predict(model, lambda m, s, rows: m(*inputs.batch(s, rows, device)), inputs.splits['val'], t['batch_size'])
        best = c.metrics(inputs.splits['val'].labels, val)
        best_state = engine.cpu_tree(model.state_dict())
        best_epoch = 0
        history = [dict(epoch=0, val_ap=best['ap'], val_auc=best['auc'], val_brier=best['brier'])]
    for epoch in range(start + 1, epochs + 1):
        clock = time.monotonic()
        model.train()
        loss_sum = 0.0
        for step, take in enumerate(c.batches(len(inputs.splits['train']), t['batch_size'], a.seed + epoch)):
            optimizer.zero_grad(set_to_none=True)
            z = model(*inputs.batch('train', take, device))
            y = torch.tensor(inputs.splits['train'].labels[take], device=device, dtype=torch.float32)
            loss = F.binary_cross_entropy_with_logits(z, y)
            require(torch.isfinite(loss), 'Nonfinite training loss')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
            optimizer.step()
            loss_sum += float(loss.detach().cpu()) * len(take)
            if step and step % 200 == 0:
                logging.info('%s epoch=%d step=%d', a.arm, epoch, step)
        val = engine.predict(model, lambda m, s, rows: m(*inputs.batch(s, rows, device)), inputs.splits['val'], t['batch_size'])
        m = c.metrics(inputs.splits['val'].labels, val)
        scheduler.step(m['ap'])
        train_z = engine.predict(model, lambda m, s, rows: m(*inputs.batch(s, rows, device)), inputs.splits['train'], t['batch_size'])
        train_eval = c.metrics(inputs.splits['train'].labels, train_z)
        train_eval_bce = float(np.mean(np.logaddexp(0, train_z) - inputs.splits['train'].labels * train_z))
        if m['ap'] > best['ap'] + 1e-06:
            best, best_epoch, best_state, stale = (m, epoch, engine.cpu_tree(model.state_dict()), 0)
        else:
            stale += 1
        elapsed += time.monotonic() - clock
        history.append(dict(epoch=epoch, train_bce=loss_sum / len(inputs.splits['train']), val_ap=m['ap'], val_auc=m['auc'], val_brier=m['brier'], train_eval_ap=train_eval['ap'], train_eval_bce=train_eval_bce, val_eval_bce=float(np.mean(np.logaddexp(0, val) - inputs.splits['val'].labels * val)), seconds=elapsed, lrs=[g['lr'] for g in optimizer.param_groups]))
        ck = dict(signature=signature, model=engine.cpu_tree(model.state_dict()), optimizer=engine.cpu_tree(optimizer.state_dict()), scheduler=scheduler.state_dict(), rng=engine.rng_state(device), epoch=epoch, stale=stale, best_state=best_state, best_metrics=best, best_epoch=best_epoch, seconds=elapsed, history=history, complete=stale >= patience or epoch == epochs)
        c.save_pt(last, ck)
        c.save_pt(best_path, dict(state=best_state, metrics=best, best_epoch=best_epoch, signature=signature))
        pd.DataFrame(history).to_csv(directory / 'history.tsv', sep='\t', index=False)
        logging.info('%s epoch=%d/%d BCE=%.5f valAP=%.6f best=%d', a.arm, epoch, epochs, history[-1]['train_bce'], m['ap'], best_epoch)
        engine.release(device)
        if ck['complete']:
            break
    model.load_state_dict(best_state)
    return (model, ck)



let plugin = LiteMol.Plugin.create({
    target: document.getElementById('litemol'),
    viewportBackground: '#ffffff',
    layoutState: { hideControls: true } // you can also include isExpanded: true
});
plugin.context.logger.message(`Hello LiteMol`);



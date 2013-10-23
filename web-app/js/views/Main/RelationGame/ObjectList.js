define(['marionette', 'templates', 'vent',
        'views/Main/RelationGame/ObjectItem'],
    function (Marionette, templates, vent,
              ObjectItem) {
  'use strict';

  return Marionette.CollectionView.extend({
    itemView: ObjectItem,
    tagName : 'ul',
    className : 'list-unstyled'

  });
});

SELECT  `ner_ann`.`text`,
        COUNT(`ner_ann`.`id`) as `occurances`

FROM `entity_recognition_entityrecognitionannotation` as `ner_ann`

INNER JOIN `document_annotation` as `ann`
  ON `ann`.`object_id` = `ner_ann`.`id` AND `ann`.`content_type_id` = (SELECT `id` FROM `django_content_type` WHERE `django_content_type`.`model` = 'entityrecognitionannotation')

INNER JOIN `document_view` as `view`
  ON `view`.`id` = `ann`.`view_id`

INNER JOIN `document_section` as `section`
  ON `section`.`id` = `view`.`section_id`

INNER JOIN `document_document` as `document`
ON `document`.`id` = `section`.`document_id` AND `document`.`id` = {document_pk}

WHERE (`ner_ann`.`type_idx` = {type_idx} AND `ner_ann`.`text` != '')

GROUP BY `ner_ann`.`text`

ORDER BY `occurances` DESC

LIMIT 20
